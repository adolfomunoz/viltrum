#pragma once
#include <random>

namespace viltrum {


//I am implementing it so it is predictable once we initialize the RNG.
template<typename S1, typename S2>
class Concat {
    S1 seq1; S2 seq2;
    
public:
    Concat(S1&& s1, S2&& s2) : seq1(std::forward<S1>(s1)), seq2(std::forward<S2>(s2)) {}
    Concat(const S1& s1, S2&& s2) : seq1(s1), seq2(std::forward<S2>(s2)) {}
    Concat(S1&& s1, const S2& s2) : seq1(std::forward<S1>(s1)), seq2(s2) {}
    Concat(const S1& s1, const S2& s2) : seq1(s1), seq2(s2) {}

    class const_iterator {
    private:
        typename S1::const_iterator i1;
        typename S1::const_iterator i1_end;
        typename S2::const_iterator i2;
        const_iterator(typename S1::const_iterator&& i1, typename S1::const_iterator&& i1e, typename S2::const_iterator&& i2) :
            i1(std::forward<typename S1::const_iterator>(i1)), i1_end(std::forward<typename S1::const_iterator>(i1e)),
            i2(std::forward<typename S2::const_iterator>(i2)) {}

        friend class Concat<S1,S2>;
    public:
        typename S1::value_type operator*() const { return (i1 != i1_end)?(*i1):(*i2); }
        const_iterator& operator++()    {  if (i1 != i1_end) ++i1; else ++i2; return (*this); }
		bool operator==(const const_iterator& that) const { return (this->i1 == that.i1) && (this->i2 == that.i2); } //Infinite list
		bool operator!=(const const_iterator& that) const { return (this->i1 != that.i1) || (this->i2 != that.i2);  } //Infinite list

    }; 

    const_iterator begin() const { return const_iterator(seq1.begin(),seq1.end(),seq2.begin()); }
    const_iterator end() const { return const_iterator(seq1.end(),seq1.end(),seq2.end());; }
};

template<typename S1,typename S2>
auto concat(S1&& seq1, S2&& seq2) {
    return Concat<std::decay_t<S1>,std::decay_t<S2>>(std::forward<S1>(seq1),std::forward<S2>(seq2));
} 

}
