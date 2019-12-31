//
// Created by amadeus on 19-11-13.
//

#ifndef PAGRAPH_SEQHELPER_HPP
#define PAGRAPH_SEQHELPER_HPP


#include <string>
#include <functional>

class SeqHelper {
private:
    static std::string testFileType(const std::string &path);
public:
    static void autoLoadFromFile(const std::string &path,
                                 std::function<void(
                                         const std::string &,
                                         const std::string &)> f);

    static void loadFromFastq(const std::string &fastqPath,
                              std::function<void(
                                      const std::string &,
                                      const std::string &,
                                      const std::string &,
                                      const std::string &)> f);

    static void loadFromFasta(const std::string &fastaPath,
                              std::function<void(
                                      const std::string &,
                                      const std::string &)> f);
};


#endif //PAGRAPH_SEQHELPER_HPP
