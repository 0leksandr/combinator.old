# combinator

```
An old code, not following standards, but interesting problem and solution.
Service, generating combinations of objects from passed collection (array/std::array/vector, etc.).
F.e.,
    OrderedCombinator([1, 2, 3], 2) = [
        [1, 2],
        [1, 3],
        [2, 3]
    ];
    ShuffledCombinator([1, 2, 3], 3) = [
        [1, 2, 3],
        [1, 3, 2],
        [2, 1, 3],
        [2, 3, 1],
        [3, 1, 2],
        [3, 2, 1]
    ];
Functions as random-access iterator, so that allows to generate arbitrary big sequenses, requiring minimum amount of memory.
Accepts collections of any objects (integers are shows as example).
Returns iterables of different types (array/std::array/vector); required type should be declared as template argument when creating combinator object.
Supports combination of pointers to original objects (if memory resourses are limited)
```
