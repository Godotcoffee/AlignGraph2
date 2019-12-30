//
// Created by amadeus on 19-5-5.
//

#ifndef PAGRAPH_MECATALIGNDATABASE_HPP
#define PAGRAPH_MECATALIGNDATABASE_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <memory>
#include "AlignInf.hpp"
#include "seq/ISeqDatabase.hpp"
#include "AbstractAlignDatabase.hpp"
#include "ParseAlignTools.hpp"

class MecatAlignDatabase : public AbstractAlignDatabase {
public:
    explicit MecatAlignDatabase(const std::string &filePath);
};


#endif //PAGRAPH_MECATALIGNDATABASE_HPP
