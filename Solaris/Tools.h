#ifndef TOOLS_H_
#define TOOLS_H_

#include <list>
#include <string>

class Component;
class Phase;

enum OS {
	Unknown,
	Windows,
	Linux
};

class Tools {
public:
	static void ToPhase(double *y, Phase *phase);
	static void	CheckAgainstSmallestNumber(const int n, double *y);

	static bool IsNumber(const std::string& s);
	static void TrimRight(std::string& str);
	static void TrimLeft(std::string& str);
	static void Trim(std::string& str);
	static void GuidToCharArray(std::string& guid, char *result);
	
	static std::string GetDirectory(const std::string path, const char directorySeparator);
	static std::string GetFileName(const std::string path, const char directorySeparator);
	static int	GetDirectorySeparator(char *c);
	static int	GetWorkingDirectory(std::string &wd);
	static enum OS GetOs();
	static int	CreatePath(const std::string &directory, const std::string &fileName, const char separator, char **path);

	static void TimeDifference(int diff, std::string& result);

	static int MergeComponentList(std::list<Component> *list1, std::list<Component> *list2, std::list<Component> *result);
	static bool Contains(std::list<Component> *list1, Component &component);
	static Component* FindByName(std::list<Component> *list, std::string name);
};

#endif
