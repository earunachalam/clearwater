#ifndef SGN_H
#define SGN_H



template <typename T>
T sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

#endif
