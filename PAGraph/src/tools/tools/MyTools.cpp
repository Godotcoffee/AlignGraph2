//
// Created by amadeus on 19-2-25.
//

#include "MyTools.hpp"
#include <iomanip>

std::chrono::milliseconds::rep MyTools::timeStamp() {
    auto d = std::chrono::steady_clock::now().time_since_epoch();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(d);
    return ms.count();
}

void MyTools::progress(double progress, int barWidth, std::ostream &of) {
    progress = std::min(1.0, progress);
    of << "\r[";
    int pos = static_cast<int>(barWidth * progress);

    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) { of << "="; }
        else if (i == pos) { of << ">"; }
        else { of << " "; }
    }

    of << "]" << std::setprecision(2) << std::fixed << (progress * 100) << "%" << " [Mem=" << totalVirtualMem() << "KB]";
    of.flush();
}
