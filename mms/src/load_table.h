#ifndef LOAD_TABLE_H
#define LOAD_TABLE_H


#include <csignal>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


#include "to_numeric.h"
#include "typedefs.h"



template <typename T>
void load_table(std::string filename, std::vector<std::vector<T> >& table, uint& nrows, uint& ncols, bool warning_ncol = true)
{
    std::ifstream file(filename);
    std::string line, token;

    nrows = 0;
    ncols = 0;

    while(std::getline(file, line))
    {
        std::istringstream iss(line);

        std::vector<T> tokens;
        while(std::getline(iss, token, ' '))
        {
            tokens.push_back(to_numeric<T>(token));
        }
        
        if (warning_ncol)
        {
            if (nrows++ == 0) ncols = tokens.size();
            else
            {
                if (ncols != tokens.size())
                {
                    fprintf(stderr, "Error: inconsistent number of columns (%zu) in file: %s.\n", tokens.size(), filename.c_str());
                    raise(SIGTRAP);
                }
            }
        }

        table.push_back(tokens);
    }
}


template <typename T>
void load_flattened_table(std::string filename, std::vector<T>& table, uint& nrows, uint& ncols, bool warning_ncol = true)
{
    std::ifstream file(filename);
    std::string line, token;

    nrows = 0;
    ncols = 0;

    while(std::getline(file, line))
    {
        std::istringstream iss(line);

        std::vector<T> tokens;
        while(std::getline(iss, token, ' '))
        {
            table.push_back(to_numeric<T>(token));
        }
        
        if (warning_ncol)
        {
            if (nrows++ == 0) ncols = tokens.size();
            else
            {
                if (ncols != tokens.size())
                {
                    fprintf(stderr, "Error: inconsistent number of columns in file.\n");
                    raise(SIGTRAP);
                }
            }
        }
    }
    std::cout << table.size() << std::endl;
}


#endif
