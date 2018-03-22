#ifndef COMBINATOR_CPP
#define COMBINATOR_CPP


#include <vector>
#include <math.h>
#include <functional>
#include <type_traits>
#include <array>
#include <boost/multi_array.hpp>
#include "my/macro.cpp"
#include "my/exception.cpp"

namespace Combinator {
    typedef unsigned long position;

    template<class element, class Container, class Combination>
    class Iterator;

    template<class element, class Container, class Combination>
    class Combinator {
            friend class Iterator<element, Container, Combination>;
        public:
            virtual Combination& operator[](position index) const = 0;
            position size() const {
                Assert(_size > 0);
                return _size;
            }
        protected:
            const Container elements;
            const position _size;
            Combinator(const Container elements, const position size):
                    elements(elements),
                    _size(size) {}
            position nrElements() const {
                return elements.size();
            }
            Combinator(const Combinator<element, Container, Combination>& other):
                    elements(other.elements),
                    _size(other._size) {}
    };

    template<class element, class Container, class Combination>
    class FixedCombinator : public Combinator<element, Container, Combination> {
            friend class Iterator<element, Container, Combination>;
        public:
            Iterator<element, Container, Combination> begin() const {
                return Iterator<element, Container, Combination>(this, 0);
            }
            Iterator<element, Container, Combination> end() const {
                return Iterator<element, Container, Combination>(this, this->size());
            }
        protected:
            const position length;
            mutable Iterator<element, Container, Combination> current;
            virtual void next(Iterator<element, Container, Combination>& iterator) const = 0;
            FixedCombinator(Container elements, const position length, const position size):
                    Combinator<element, Container, Combination>(elements, size),
                    length(length),
                    current(this, 0) {}
            void first(position* const positions) const {
                for (position c = 0; c < this->length; c++)
                    positions[c] = c;
            }
            FixedCombinator(const FixedCombinator<element, Container, Combination>& other):
                    Combinator<element, Container, Combination>(other),
                    length(other.length),
                    current(this, 0) {}
    };

    template<class element, class Container, class Combination = Container>
    class OrderedCombinator;
    template<class element, class Container, class Combination = Container>
    class ShuffledCombinator;

    template<class element, class Container, class Combination>
    class Iterator {
            friend class OrderedCombinator<element, Container, Combination>;
            friend class ShuffledCombinator<element, Container, Combination>;
            class Converter;
        public:
            Iterator(
                    const FixedCombinator<element, Container, Combination>* const combinator,
                    const position index
            ):
                    combinator(combinator),
                    index(index),
                    positions(new position[combinator->length]),
                    converter(combinator),
                    combination(converter.construct(&combination)) {
                combinator->first(positions);
            }
            ~Iterator() {
                delete positions;
                converter.destruct(combination);
            }
            void operator++() {
                combinator->next(*this);
                index++;
            }
            bool operator!=(const Iterator<element, Container, Combination>& other) const {
                return index != other.index;
            }
            Combination& operator*() {
                converter.prepare(combination);
                for (position c = 0; c < combinator->length; c++)
                    combination[c] = converter.getElement(
                            &combination,
                            combinator->elements[positions[c]]
                    );
                return combination;
            }
        protected:
            position index;
            position* const positions;
        private:
            const FixedCombinator<element, Container, Combination>* const combinator;
            const Converter converter;
            Combination combination;

            class Converter {
                public:
                    Converter(const FixedCombinator<
                            element,
                            Container,
                            Combination
                    >* const combinator):
                            combinator(combinator) {}
                    std::vector<element> construct(std::vector<element>*) const {
                        std::vector<element> vec;
                        initVector(vec);
                        return vec;
                    }
                    template<class _element>
                    simpleArray<_element> construct(simpleArray<_element>*) const {
                        return simpleArray<_element>(combinator->length);
                    }
                    template<class _element, unsigned long size>
                    std::array<_element, size> construct(std::array<_element, size>*) const {
                        Assert(size == combinator->length);
                        return std::array<_element, size>();
                    }
                    template<class _element>
                    _element* construct(_element**) const {
                        return new _element[combinator->length];
                    }
                    template<class C>
                    C construct(C*) const {
                        throw exception("Unsupported combination type");
                    }

                    void prepare(std::vector<element>& combination) const {
                        if (combination.size() != combinator->length) {
                            combination.clear();
                            initVector(combination);
                        }
                    }
                    template<class C>
                    void prepare(C& combination) const {}

                    element* getElement(simpleArray<element*>*, element& _element) const {
                        return &_element;
                    }
                    template<unsigned long size>
                    element* getElement(std::array<element*, size>*, element& _element) const {
                        return &_element;
                    }
                    element* getElement(element***, element& _element) const {
                        return &_element;
                    }
                    template<class C>
                    const element& getElement(C*, const element& _element) const {
                        return _element;
                    }

                    template<class _element>
                    void destruct(_element*& combination) const {
                        delete combination;
                    }
                    template<class C>
                    void destruct(C& combination) const {}
                private:
                    const FixedCombinator<element, Container, Combination>* const combinator;
                    void initVector(std::vector<element>& vec) const {
                        vec.reserve(combinator->length);
                        for (position c = 0; c < combinator->length; c++)
                            vec.push_back(combinator->elements[c]);
                    }
            };
    };

    template<class element, class Container, class Combination>
    class OrderedCombinator : public FixedCombinator<element, Container, Combination> {
            class OrderIterator;
            class Hunter;
            class Mathematician;
        public:
            OrderedCombinator(Container elements, const position length):
                    FixedCombinator<element, Container, Combination>(
                            elements,
                            length,
                            this->nPerM(elements.size(), length)
                    ),
                    iterators() {}
            ~OrderedCombinator() {
                for (OrderIterator* iterator : iterators)
                    delete iterator;
            }
            Combination& operator[](position index) const override {
                if (iterators.size() == 0)
                    initIterators();
                position estimated((position)-1);
                OrderIterator* chosen(NULL);
                for (OrderIterator* iterator : iterators) {
                    position myBet = iterator->estimate(index);
                    if (myBet < estimated) {
                        chosen = iterator;
                        estimated = myBet;
                    }
                }
                chosen->go(index);
                for (position c = 0; c < this->length; c++)
                    this->current.positions[c] = chosen->positions[c];
                return *this->current;
            }
            OrderedCombinator(const OrderedCombinator<element, Container, Combination>& other):
                    FixedCombinator<element, Container, Combination>(other),
                    iterators() {}
        protected:
            void next(Iterator<element, Container, Combination>& iterator) const override {
                next(iterator.positions);
            }
        private:
            mutable std::vector<OrderIterator*> iterators;
            void initIterators() const {
                iterators.push_back(new Hunter(this));
                iterators.push_back(new Mathematician(this));
            }
            void increment(position* const positions, const position _position) const {
                if (++positions[_position] > maxPosition(_position)) {
                    if (_position != 0) { // preventing fail on end
                        increment(positions, _position - 1);
                        positions[_position] = positions[_position - 1] + 1;
                    }
                }
            }
            void decrement(position* const positions, const position _position) const {
                --positions[_position];
                if (_position > 0 && positions[_position] == positions[_position - 1]) {
                    for (position c = _position; c < this->length; c++)
                        positions[c] = maxPosition(c);
                    decrement(positions, _position - 1);
                }
            }
            void next(position* const positions) const {
                increment(positions, this->length - 1);
            }
            void previous(position* const positions) const {
                decrement(positions, this->length - 1);
            }
            position maxPosition(const position _position) const { // TODO: inline?
                return this->nrElements() + _position - this->length;
            }
            position nPerM(const position n, const position m) const {
                double res(1.);
                for (position c = 0; c < m; c++) {
                    res *= n - c;
                    res /= c + 1;
                }
                return (position)res;
            }

            class OrderIterator {
                public:
                    OrderIterator(const OrderedCombinator<
                            element,
                            Container,
                            Combination
                    >* const combinator):
                            positions(new position[combinator->length]),
                            combinator(combinator) {}
                    ~OrderIterator() {
                        delete positions;
                    }
                    position* const positions;
                    virtual position estimate(position index) const = 0;
                    virtual void go(position index) = 0;
                protected:
                    const OrderedCombinator<element, Container, Combination>* const combinator;
            };
            class Walker : public OrderIterator {
                public:
                    Walker(const OrderedCombinator<
                            element,
                            Container,
                            Combination
                    >* const combinator):
                            OrderIterator(combinator),
                            location(0) {
                        this->combinator->first(this->positions);
                    }
                    position location;
                    position estimate(const position index) const override {
                        if (index > location)
                            return index - location;
                        else if (index < location)
                            return location - index;
                        else
                            return 0;
                    }
                    void go(const position index) override {
                        while (location < index)
                            forward();
                        while (location > index)
                            back();
                    }
                    Walker(const Walker& other):
                            OrderIterator(other.combinator),
                            location(other.location) {
                        for (position c = 0; c < this->combinator->length; c++)
                            this->positions[c] = other.positions[c];
                    }
                    void forward() {
                        this->combinator->next(this->positions);
                        ++location;
                    }
                    void back() {
                        this->combinator->previous(this->positions);
                        --location;
                    }
            };
            class Hunter : public OrderIterator {
                public:
                    Hunter(const OrderedCombinator<
                            element,
                            Container,
                            Combination
                    >* const combinator):
                            OrderIterator(combinator),
                            nrElements(combinator->nrElements()) {
                        const position nrGuardians = (position)sqrt(this->combinator->size()) + 1;
                        reactionTime = this->combinator->size() / nrGuardians; // TODO: check
                        Walker patrol(this->combinator);
                        while (patrol.location < this->combinator->size() - 1) {
                            patrol.forward();
                            if ((patrol.location + reactionTime / 2) % reactionTime == 0)
                                guardians.push_back(patrol);
                        }
                    }
                    position estimate(const position index) const override {
                        return guardian(index).estimate(index);
                    }
                    void go(const position index) override {
                        Walker envoy = guardian(index);
                        envoy.go(index);
                        for (position c = 0; c < this->combinator->length; c++)
                            this->positions[c] = envoy.positions[c];
                    }
                private:
                    std::vector<Walker> guardians; // TODO: use array?
                    const position nrElements;
                    position reactionTime;
                    Walker guardian(const position index) const {
                        return guardians[index / reactionTime];
                    }
            };
            class Mathematician : public OrderIterator {
                public:
                    Mathematician(const OrderedCombinator<
                            element,
                            Container,
                            Combination
                    >* const combinator):
                            OrderIterator(combinator) {
                        avgEstimation = avgNrSteps();
                    }
                    position estimate(position index) const override {
                        return avgEstimation; // TODO: test
                    }
                    void go(position index) override {
                        position nrElements(this->combinator->nrElements());
                        for (position c = 0; c < this->combinator->length; c++) {
                            step _step = getStep(c, nrElements, index);
                            nrElements -= _step.x + 1;
                            index -= _step.beginningOfX;
                            insertUnique(c, _step.x);
                        }
                    }
                private:
                    position avgEstimation;
                    struct step {
                        step(
                                const position x,
                                const position beggingOfX
                        ): x(x), beginningOfX(beggingOfX) {}
                        position x;
                        position beginningOfX;
                    };
                    step getStep(
                            const position _position,
                            const position nrElements,
                            const position index
                    ) const {
                        if (_position == this->combinator->length - 1)
                            return step(index, index);
                        step res(0, 0);
                        position lastNPerM(0);
                        while (res.beginningOfX <= index) { // TODO: optimize
                            res.beginningOfX += lastNPerM = this->combinator->nPerM(
                                    nrElements - res.x - 1,
                                    this->combinator->length - _position - 1
                            );
                            res.x++;
                        }
                        if (res.x > 0)
                            --res.x;
                        res.beginningOfX -= lastNPerM;
                        return res;
                    }
                    void insertUnique(const position _position, position value) const {
                        if (_position > 0)
                            value += this->positions[_position - 1] + 1;
                        this->positions[_position] = value;
                    }
                    position avgNrSteps() const {
                        float totalAvgNrSteps(0.f);
                        position prevValue(0);
                        for (position c = 0; c < this->combinator->length; c++) {
                            float nrSteps = avgNrSteps(c, prevValue);
                            prevValue = getStep(
                                    c,
                                    this->combinator->nrElements() - prevValue,
                                    (position)nrSteps
                            ).x;
                            totalAvgNrSteps += nrSteps;
                        }
                        return (position)totalAvgNrSteps;
                    }
                    float avgNrSteps(
                            const position _position,
                            const position minValue
                    ) const {
                        float sum(0.f);
                        position maxValue = this->combinator->nrElements()
                                            + _position
                                            + 1
                                            - this->combinator->length;
                        for (position c = minValue; c < maxValue; c++)
                            sum += this->combinator->nPerM(
                                    this->combinator->nrElements() - c - 1,
                                    this->combinator->length - _position - 1
                            );
                        return sum / (maxValue - minValue + 1);
                    }
            };
    };

    template<class element, class Container, class Combination>
    class ShuffledCombinator : public FixedCombinator<element, Container, Combination> {
        public:
            ShuffledCombinator(Container elements, const position length):
                    FixedCombinator<element, Container, Combination>(
                            elements,
                            length,
                            this->nPerM(elements.size(), length)
                    ) {}
            Combination& operator[](const position index) const override {
                move(this->current, index);
                return *this->current;
            }
        protected:
            void next(Iterator<element, Container, Combination>& iterator) const override {
                move(iterator, iterator.index + 1);
            }
        private:
            void move(
                    const Iterator<element, Container, Combination>& iterator,
                    position index
            ) const {
                position nrElements(this->nrElements());
                for (position c = 0; c < this->length; c++) {
                    insertUnique(iterator, c, index % nrElements);
                    index /= nrElements--;
                }
            }
            void insertUnique(
                    const Iterator<element, Container, Combination>& iterator,
                    const position _position,
                    position value
            ) const {
                position prevValue((position)-1);
                do {
                    position add(0); // TODO: move up?
                    for (position c = 0; c < _position; c++)
                        if (
                                prevValue + 1 < iterator.positions[c] + 1
                                && iterator.positions[c] <= value)
                            ++add;
                    prevValue = value;
                    value += add;
                } while (prevValue != value);
                iterator.positions[_position] = value;
            }
            position nPerM(const position n, const position m) const {
                return nPerM0(n, m);
            }
            position nPerM0(const position n, const position m) const {
                if (m > 1)
                    return n * nPerM0(n - 1, m - 1);
                return n;
            }
            position nPerM1(position n, position m) const {
                position res(n);
                while (m --> 1)
                    res *= --n;
                return res;
            }
            position nPerM2(position n, position m) const {
                position res(n);
                while (--m > 0)
                    res *= --n;
                return res;
            }
    };

    namespace tests {
        void testOrdered() {
            std::vector<double> vec;
            for (unsigned c = 0; c < 4; c++)
                vec.push_back(c + 1);

            const unsigned NR_ELEMENTS_IN_COMBINATION = 2;
            OrderedCombinator<double, std::vector<double>, std::vector<double>> combinations(
                    vec,
                    NR_ELEMENTS_IN_COMBINATION
            );
            const unsigned NR_COMBINATIONS = 6;
            double expectedOrdered[NR_COMBINATIONS][2] = {
                    {1., 2.},
                    {1., 3.},
                    {1., 4.},
                    {2., 3.},
                    {2., 4.},
                    {3., 4.},
            };
            Assert(combinations.size() == NR_COMBINATIONS);
            position c(0);
            for (auto combination : combinations) {
                Assert(combination.size() == NR_ELEMENTS_IN_COMBINATION);
                Assert(combinations[c].size() == NR_ELEMENTS_IN_COMBINATION);
                for (unsigned d = 0; d < NR_ELEMENTS_IN_COMBINATION; d++) {
                    double expected = expectedOrdered[c][d];
                    Assert(combination[d] == expected);
                    Assert(combinations[c][d] == expected);
                }
                c++;
            }
            myPrint("Test passed\n");
        }
        void testShuffled() {
            std::vector<double> vec;
            for (unsigned c = 0; c < 4; c++)
                vec.push_back(c + 1);

            const unsigned NR_ELEMENTS_IN_COMBINATION = 2;
            ShuffledCombinator<double, std::vector<double>, std::vector<double>> combinations(
                    vec,
                    NR_ELEMENTS_IN_COMBINATION
            );
            const unsigned NR_COMBINATIONS = 12;
            double expectedShuffled[NR_COMBINATIONS][2] = {
                    {1., 2.},
                    {2., 1.},
                    {3., 1.},
                    {4., 1.},
                    {1., 3.},
                    {2., 3.},
                    {3., 2.},
                    {4., 2.},
                    {1., 4.},
                    {2., 4.},
                    {3., 4.},
                    {4., 3.},
            };
            Assert(combinations.size() == NR_COMBINATIONS);
            position c(0);
            for (auto combination : combinations) {
                Assert(combination.size() == NR_ELEMENTS_IN_COMBINATION);
                Assert(combinations[c].size() == NR_ELEMENTS_IN_COMBINATION);
                for (unsigned d = 0; d < NR_ELEMENTS_IN_COMBINATION; d++) {
                    double expected = expectedShuffled[c][d];
                    Assert(combination[d] == expected);
                    Assert(combinations[c][d] == expected);
                }
                c++;
            }
            myPrint("Test passed\n");
        }
        void tests() {
            testOrdered();
            testShuffled();
        }
    }
}


#endif //COMBINATOR_CPP
