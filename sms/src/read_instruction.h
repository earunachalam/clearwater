#ifndef READ_INSTRUCTION_H
#define READ_INSTRUCTION_H


#include <string>
#include <vector>



template <typename T>
void read_instruction(std::vector<std::string>&);

template <typename T>
void read_instruction(std::vector<std::string>&, std::vector<T*>&);


#endif
