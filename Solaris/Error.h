#ifndef ERROR_H_
#define ERROR_H_

#include <list>
#include <string>

class Error {
public:
	static void PushLocation(const char *file, const char *function, long line);
	static void PrintStackTrace();

	static std::string _errMsg;

private:
	static std::list<std::string> _stackTrace;
};

#endif
