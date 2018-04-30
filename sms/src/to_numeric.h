#ifndef TO_NUMERIC_H
#define TO_NUMERIC_H


#include <csignal>
#include <string>


// convert string to numeric type
// primary template and specializations to different numeric datatypes


template <typename T>
T to_numeric(std::string s)
{
    return T();
}


template <>
unsigned int to_numeric<unsigned int>(std::string s)
{
    return static_cast<unsigned int>(std::stoi(s));
}

template <>
int to_numeric<int>(std::string s)
{
    return std::stoi(s);
}

template <>
double to_numeric<double>(std::string s)
{
    return std::stod(s);
}



#endif
