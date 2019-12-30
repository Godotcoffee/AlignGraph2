//
// Created by amadeus on 19-6-24.
//

#ifndef PAGRAPH_MUMMERALIGNDATABASEV2_HPP
#define PAGRAPH_MUMMERALIGNDATABASEV2_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <memory>
#include "AlignInf.hpp"
#include "AbstractAlignDatabase.hpp"
#include "ParseAlignTools.hpp"

class MummerAlignDatabaseV2 : public AbstractAlignDatabase {
public:
    explicit MummerAlignDatabaseV2(const std::string &filePath);
};


#endif //PAGRAPH_MUMMERALIGNDATABASEV2_HPP
