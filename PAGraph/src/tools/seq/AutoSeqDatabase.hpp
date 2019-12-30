//
// Created by amadeus on 19-12-26.
//

#ifndef PAGRAPH_AUTOSEQDATABASE_HPP
#define PAGRAPH_AUTOSEQDATABASE_HPP


#include "AbstractSeqDatabase.hpp"

class AutoSeqDatabase : public AbstractSeqDatabase {
public:
    explicit AutoSeqDatabase(const std::string &filePath);
};


#endif //PAGRAPH_AUTOSEQDATABASE_HPP
