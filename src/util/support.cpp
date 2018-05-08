/** 
 * @file support.cpp
 * @author Arne Kutzner
 */
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

#include "util/support.h"
#include "util/exception.h"
/* Constructs the full file name for some prefix, suffix combination.
 * Returns by value for convenience purposes.
 */
std::string fullFileName( const char *pcFileNamePrefix, const char *pcSuffix )
{
    std::string sFileName( pcFileNamePrefix );
    sFileName.push_back('.');
    return sFileName.append( pcSuffix );
} // method

bool fileExists(const std::string& rsFile)
{
    struct stat buffer;
    return (stat (rsFile.c_str(), &buffer) == 0);
}// function

void makeDir(const std::string& rsFile)
{
    mode_t nMode = 0733; // UNIX style permissions
    int nError = 0;
    #if defined(_WIN32)
        nError = _mkdir(rsFile.c_str()); // can be used on Windows
    #else 
        nError = mkdir(rsFile.c_str(), nMode); // can be used on non-Windows
    #endif
    if (nError != 0) {
        throw NullPointerException(std::string("Could not create Dir: ").append(rsFile).c_str());
    }// if
}// function
