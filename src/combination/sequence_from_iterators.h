#pragma once

namespace viltrum {

template<typename Iterator>
class sequence_from_iterators {
    Iterator ibegin, iend;
public:
    using const_iterator = Iterator;
    
    sequence_from_iterators(const Iterator& ibegin, const Iterator& iend) : ibegin(ibegin), iend(iend) {}

    Iterator begin() const { return ibegin; }
    Iterator end() const { return end; }
};

}