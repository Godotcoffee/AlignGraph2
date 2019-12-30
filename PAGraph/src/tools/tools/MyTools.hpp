//
// Created by amadeus on 19-2-25.
//

#ifndef PAGRAPH_MYTOOLS_HPP
#define PAGRAPH_MYTOOLS_HPP

#include <iostream>
#include <cmath>
#include <chrono>
#include <memory>
#include <cstring>
#include <cstdlib>


#include "sys/types.h"
#include "sys/sysinfo.h"

class MyTools {
private:
    static long long parseLine(char* line){
        // This assumes that a digit will be found and the line ends in " Kb".
        long long  i = strlen(line);
        const char* p = line;
        while (*p <'0' || *p > '9') p++;
        line[i-3] = '\0';
        i = atoll(p);
        return i;
    }

    static long long getValue(){ //Note: this value is in KB!
        FILE* file = fopen("/proc/self/status", "r");
        long long result = -1;
        char line[128];

        while (fgets(line, 128, file) != nullptr){
            if (strncmp(line, "VmSize:", 7) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    }
public:
    static std::chrono::milliseconds::rep timeStamp();

    static void progress(double progress, int barWidth, std::ostream &of);

    template<typename T, typename... Args>
    static std::unique_ptr<T> make_unique(Args&&... args)
    {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    static long long totalVirtualMem() {
        return getValue();
    }


};


#endif //PAGRAPH_MYTOOLS_HPP
