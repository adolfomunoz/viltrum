/**
 * Developed within the framework of Diego Viñuales' final bachelor's degree project at the University of Zaragoza,
 * under the supervision of Dr. Adolfo Muñoz
 * 
 * From repository:
 * Relaxed Priority Queues - https://github.com/DiegoV8/relaxedPriorityQueues_874280.git
 */

#ifndef SKIPLIST_HPP
#define SKIPLIST_HPP

#include <vector>
#include <optional>
#include <functional>
#include <mutex>
#include <atomic>
#include <thread>
#include <array>
#include <random>
#include <algorithm>
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
#include <immintrin.h>
#endif

// Clase Spinlock usada para gestionar el bloqueo de nodos de Skiplist.
class Spinlock {
    std::atomic_flag locked = ATOMIC_FLAG_INIT;
public:
    void lock() {
        int backoff = 1;
        while (locked.test_and_set(std::memory_order_acquire)) {
            for (int i = 0; i < backoff; ++i) {
                #if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
                                _mm_pause();
                #else
                                std::this_thread::yield();
                #endif
            }
            if (backoff < 1024) backoff *= 2;
            else std::this_thread::yield();
        }
    }

    bool try_lock() {
        return !locked.test_and_set(std::memory_order_acquire);
    }

    void unlock() {
        locked.clear(std::memory_order_release);
    }
};

/**
 * @brief Interfaz para una Skiplist.
 * @tparam T Tipo de dato.
 * @tparam Compare Comparador.
 */
template <typename T, typename Compare = std::greater<T>>
class Skiplist {
public:
    using value_type = T;

    /**
     * @brief Constructor para la Skiplist.
     * @param comp Instancia del comparador de prioridad.
     * @param rel_factor Grado de relajación de la estructura.
     */
    explicit Skiplist(Compare comp = Compare(), int rel_factor = 0)
        : current_level(1), comp(comp), relaxation_factor(rel_factor)
    {
        NIL    = new node(std::size_t(0), nullptr);
        header = new node(MAX_LEVEL, NIL);
    }

    /**
     * @brief Destructor para la Skiplist.
     */
    ~Skiplist() {
        node* curr = header->forward[0].load(std::memory_order_relaxed);
        while (curr != NIL) {
            node* next = curr->forward[0].load(std::memory_order_relaxed);
            if (!curr->marked_for_deletion.load(std::memory_order_relaxed)) {
                curr->get_data().~T();
                delete curr;
            }
            curr = next;
        }
        delete header;
        delete NIL;

        for (node* n : retired_nodes)
            delete n;
    }

    /**
     * @brief Inserta un elemento de forma concurrente.
     * @param value Valor a insertar en la estructura
     */
    void push(T value) {
        std::size_t new_level = random_level();
        node* new_node = new node(value, new_level, NIL);

        std::array<node*, MAX_LEVEL> preds;
        std::array<node*, MAX_LEVEL> succs;

        while (true) {
            node* curr = header;
            for (int i = MAX_LEVEL - 1; i >= 0; --i) {
                node* next = curr->forward[i].load(std::memory_order_acquire);
                while (next != NIL &&
                       !next->marked_for_deletion.load(std::memory_order_acquire) &&
                       comp(next->get_data(), value)) {
                    curr = next;
                    next = curr->forward[i].load(std::memory_order_acquire);
                }
                preds[i] = curr;
                succs[i] = next;
            }

            std::vector<node*> to_lock;
            to_lock.reserve(new_level);
            for (std::size_t i = 0; i < new_level; ++i) {
                if (std::find(to_lock.begin(), to_lock.end(), preds[i]) == to_lock.end())
                    to_lock.push_back(preds[i]);
            }
            std::sort(to_lock.begin(), to_lock.end());
            for (node* n : to_lock) n->node_lock.lock();

            bool failed = false;
            for (std::size_t i = 0; i < new_level; ++i) {
                if (preds[i]->marked_for_deletion.load(std::memory_order_acquire) ||
                    preds[i]->forward[i].load(std::memory_order_acquire) != succs[i]) {
                    failed = true;
                    break;
                }
            }

            if (failed) {
                for (node* n : to_lock) n->node_lock.unlock();
                continue;
            }

            for (std::size_t i = 0; i < new_level; ++i) {
                new_node->forward[i].store(succs[i], std::memory_order_relaxed);
                preds[i]->forward[i].store(new_node, std::memory_order_release);
            }

            new_node->fully_linked.store(true, std::memory_order_release);
            for (node* n : to_lock) n->node_lock.unlock();
            break;
        }

        _size.fetch_add(1, std::memory_order_relaxed);
    }

    /**
     * @brief Intenta extraer el elemento de mayor prioridad.
     * @return std::optional con el valor, o nullopt si la estructura está saturada o vacía.
     */
    std::optional<T> try_pop() {
        node* target = nullptr;
        T extracted_value;

        std::size_t current_size = _size.load(std::memory_order_relaxed);
        int steps = (relaxation_factor > 0 && current_size > (std::size_t)relaxation_factor)
                    ? (int)fast_rand(relaxation_factor) : 0;

        node* curr = header->forward[0].load(std::memory_order_acquire);
        while (curr != NIL) {
            if (curr->fully_linked.load(std::memory_order_acquire) &&
                !curr->marked_for_deletion.load(std::memory_order_acquire)) {

                if (steps <= 0) {
                    bool expected = false;
                    if (curr->marked_for_deletion.compare_exchange_strong(expected, true, std::memory_order_acquire)) {
                        extracted_value = curr->get_data();
                        target = curr;
                        break;
                    }
                } else {
                    steps--;
                }
            }
            curr = curr->forward[0].load(std::memory_order_acquire);
        }

        if (target == nullptr) {
            return std::nullopt;
        }

        for (int i = MAX_LEVEL - 1; i >= 0; --i) {
            bool unlinked = false;
            
            while (!unlinked) {
                node* pred = header;
                node* next = pred->forward[i].load(std::memory_order_acquire);

                while (next != NIL && next != target) {
                    if (comp(next->get_data(), extracted_value) || !comp(extracted_value, next->get_data())) {
                        if (!next->marked_for_deletion.load(std::memory_order_acquire)) {
                            pred = next;
                        }
                        next = next->forward[i].load(std::memory_order_acquire);
                    } else {
                        break;
                    }
                }

                if (next == target) {
                    pred->node_lock.lock();
                    if (!pred->marked_for_deletion.load(std::memory_order_relaxed)) {
                        bool valid = true;
                        node* curr_trace = pred->forward[i].load(std::memory_order_relaxed);
                        
                        while (curr_trace != target) {
                            if (curr_trace == NIL || !curr_trace->marked_for_deletion.load(std::memory_order_relaxed)) {
                                valid = false;
                                break;
                            }
                            curr_trace = curr_trace->forward[i].load(std::memory_order_relaxed);
                        }
                        
                        if (valid) {
                            node* target_next = target->forward[i].load(std::memory_order_relaxed);
                            pred->forward[i].store(target_next, std::memory_order_release);
                            unlinked = true;
                        }
                    }
                    pred->node_lock.unlock();
                } else {
                    unlinked = true; 
                }
            }
        }

        _size.fetch_sub(1, std::memory_order_relaxed);
        target->get_data().~T();
        retire_node(target);

        return extracted_value;
    }

    /**
     * @brief Comprueba si la estructura está vacía.
     * @return bool con True si la estructura esta vacía y False si no.
     */
    bool empty() const {
        return _size.load(std::memory_order_relaxed) == 0;
    }

    /**
     * @brief Retorna el número aproximado de elementos.
     * @return std::size_t con el valor aproximado del numero de elementos en la Skiplist.
     */
    std::size_t size() const {
        return _size.load(std::memory_order_relaxed);
    }

    /**
     * @brief Devuelve la concatenacion de todos los nodos para uso más eficiente en regions-adaptative-multiqueue. No modifica la cola.
     * @return std::vector<T> con todos los elementos en la estructura sin ordenar.
     */
    std::vector<T> drain() {
        std::vector<T> result;

        node* curr = header->forward[0].load(std::memory_order_relaxed);
        while (curr != NIL) {
            node* next = curr->forward[0].load(std::memory_order_relaxed);
            if (!curr->marked_for_deletion.load(std::memory_order_relaxed) &&
                 curr->fully_linked.load(std::memory_order_relaxed)) {
                result.push_back(curr->get_data());
                curr->marked_for_deletion.store(true, std::memory_order_relaxed);
                _size.fetch_sub(1, std::memory_order_relaxed);
                curr->get_data().~T();
                retire_node(curr);
            }
            curr = next;
        }

        return result;
    }

    private:
    static constexpr std::size_t MAX_LEVEL = 16; 

    // Estructura usada para representar cada nodo en la Skiplist.
    struct node {
        alignas(T) unsigned char data_buffer[sizeof(T)];

        std::array<std::atomic<node*>, MAX_LEVEL> forward;
        Spinlock node_lock;
        std::atomic<bool> fully_linked{false};
        std::atomic<bool> marked_for_deletion{false};

        node(T val, std::size_t /*level*/, node* nil_ptr) {
            new (data_buffer) T(std::move(val));
            for (std::size_t i = 0; i < MAX_LEVEL; ++i)
                forward[i].store(nil_ptr, std::memory_order_relaxed);
        }

        node(std::size_t /*level*/, node* nil_ptr) {
            for (std::size_t i = 0; i < MAX_LEVEL; ++i)
                forward[i].store(nil_ptr, std::memory_order_relaxed);
        }

        T& get_data() {
            return *reinterpret_cast<T*>(data_buffer);
        }
    };

    /**
    * @brief Genera un numero aleatorio en [0, max_val] usando un algoritmo Xorshift.
    * @param max_val Valor máximo generado.
    * @return std::size_t con valor en [0, max_val].
    */
    inline std::size_t fast_rand(std::size_t max_val) {
        static thread_local uint32_t state = []() {
            uint32_t seed = std::random_device{}();
            return seed == 0 ? 1 : seed;
        }();
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state % max_val;
    }

    // Atributos internos de la Skiplist.
    node* header;
    node* NIL;
    std::size_t              current_level;
    Compare                  comp;
    std::atomic<std::size_t> _size{0};
    int                      relaxation_factor;
    static constexpr float P = 0.5f;

    std::mutex           retired_mutex;
    std::vector<node*>   retired_nodes;

    /**
    * @brief Genera un nivel aleatorio para un nodo, hay una probabilidad 1/(2i−1) de alcanzar el nivel i.
    * @return std::size_t con valor en [1, MAX_LEVEL].
    */
    std::size_t random_level() {
        std::size_t lvl = 1;
        while (fast_rand(100) < (std::size_t)(P * 100) && lvl < MAX_LEVEL)
            lvl++;
        return lvl;
    }

    /**
    * @brief Manda un nodo a la lista de nodos a ser borrados.
    * @param n nodo a ser retirado.
    */
    void retire_node(node* n) {
        std::lock_guard<std::mutex> lk(retired_mutex);
        retired_nodes.push_back(n);
    }
};

#endif