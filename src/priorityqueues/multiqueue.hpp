/**
 * Developed within the framework of Diego Viñuales' final bachelor's degree project at the University of Zaragoza,
 * under the supervision of Dr. Adolfo Muñoz
 * 
 * The only change is that the "drain" method 
 * is destructive and moves elements instead of copying them. Still 
 * the improvement on speed has been small
 * 
 * From repository:
 * Relaxed Priority Queues - https://github.com/DiegoV8/relaxedPriorityQueues_874280.git
 */

#ifndef MULTI_QUEUE_HPP
#define MULTI_QUEUE_HPP

#include <vector>
#include <optional>
#include <functional>
#include <random>
#include <mutex>
#include <memory>
#include <algorithm>
#include <iterator>
#include <queue>

/**
 * @brief Interfaz para una Multiqueue.
 * @tparam T Tipo de dato.
 * @tparam Compare Comparador.
 */
template <typename T, typename Compare = std::less<T>>
class Multiqueue {
public:
    using value_type = T;

    /**
     * @brief Constructor para la Multiqueue.
     * @param n Número de sub-colas internas.
     * @param c Número de sub-colas que se comparan al extraer elementos.
     * @param comp Instancia del comparador de prioridad.
     */
    explicit Multiqueue(std::size_t n = 8, int c = 2, Compare comp = Compare()) {
        this->n = n;
        this->c = c;
        this->comp = comp;

        queues.reserve(n);

        for (std::size_t i = 0; i < n; ++i) {
            queues.push_back(std::make_unique<Q>(comp));
        }
    }

    /**
     * @brief Inserta un elemento de forma concurrente.
     * @param value Valor a insertar en la estructura
     */
    void push(const T& value){
        
        std::size_t i = fast_rand(n);
        {
            std::lock_guard<std::mutex> lock(queues[i]->M); 
            queues[i]->pq.push(value);
        }
    }

    /**
     * @brief Intenta extraer el elemento de mayor prioridad.
     * @return std::optional con el valor, o nullopt si la estructura está saturada o vacía.
     */
    std::optional<T> try_pop() {

        std::vector<std::size_t> locked_indices;
        int best_idx = -1;

        for (int i = 0; i < c; ++i) {
            std::size_t idx = fast_rand(n);
            
            if (std::find(locked_indices.begin(), locked_indices.end(), idx) != locked_indices.end()) {
                continue;
            }
            
            if (queues[idx]->M.try_lock()) {
                if (!queues[idx]->pq.empty()) {
                    locked_indices.push_back(idx);
                    
                    if (best_idx == -1 || comp(queues[best_idx]->pq.top(), queues[idx]->pq.top())) {
                        best_idx = idx;
                    }
                } else {
                    queues[idx]->M.unlock();
                }
            }
        }

        if (best_idx == -1) {
            for (std::size_t idx : locked_indices) queues[idx]->M.unlock();
            return std::nullopt;
        }

        T result = queues[best_idx]->pq.top();
        queues[best_idx]->pq.pop();

        for (std::size_t idx : locked_indices) {
            queues[idx]->M.unlock();
        }
        return result;
    }

    /**
     * @brief Comprueba si todas las colas internas están vacías.
     * @return bool con True si la estructura esta vacía y False si no.
     */
    bool empty() const {
        for (std::size_t i = 0; i < n; i++){
            std::lock_guard<std::mutex> lock(queues[i]->M); 
            if (!queues[i]->pq.empty()) return false;
        }
        return true;
    }

    /**
     * @brief Retorna el número aproximado de elementos en la estructura.
     * @return std::size_t con el valor aproximado del numero de elementos en la Multiqueue.
     */
    std::size_t size() const {
        std::size_t sum = 0;
        for (std::size_t i = 0; i < n; i++){
            std::lock_guard<std::mutex> lock(queues[i]->M); 
            sum += queues[i]->pq.size();
        }
        return sum;
    }

    /**
     * @brief Mueve todos los elementos a un vector, dejando la cola vacía. No es thread-safe.
     * @return std::vector<T> con todos los elementos en la estructura sin ordenar.
     */
    std::vector<T> drain() {
        std::vector<T> res;
        res.reserve(size());

        for (auto& q : queues) {
            auto& v = q->pq.container();
            res.insert(res.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
            v.clear();
        }
        return res;
    }
private:

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

    // Estructura usada para representar cada una de las colas en la Multiqueue.
    struct Q {
        // Expone el contenedor protegido de std::priority_queue.
        struct Heap : std::priority_queue<T, std::vector<T>, Compare> {
            using Base = std::priority_queue<T, std::vector<T>, Compare>;
            using Base::Base;
            std::vector<T>& container() { return this->c; }
        };
        Heap pq;
        mutable std::mutex M;
        Q(const Compare& comp) : pq(comp) {}
    };

    // Atributos internos de la Multiqueue.
    std::vector<std::unique_ptr<Q>> queues;
    std::size_t n;
    int c;
    Compare comp;
    
};

#endif
