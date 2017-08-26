#ifndef ARGUMENTHANDLER_H
#define ARGUMENTHANDLER_H

#include <string>
#include <vector>
#include <map>

enum ArgumentFlag {REQUIRED_ARGUMENT, OPTIONAL_ARGUMENT, SWITCH, FLEXIBLE_SWITCH};

class ArgumentWrapperInterface {
public:
	virtual ~ArgumentWrapperInterface() {}
	virtual void SetToDefaultValue() = 0;
	virtual void SetToSwitchValue() = 0;
	virtual void SetToValueFromString(std::string value) = 0;
	virtual std::string GetDescription() = 0;
	virtual ArgumentFlag GetFlag() = 0;
};

template <typename T>
class ArgumentWrapper : public ArgumentWrapperInterface {
public:
	ArgumentWrapper(T& handle) : handle(&handle) { }
	void SetValue(T value) { *handle = value; }
	ArgumentFlag flag;
	T defaultValue;
	T switchValue;
	std::string description;
	
	void SetToDefaultValue();
	void SetToSwitchValue();
	void SetToValueFromString(std::string value);
	std::string GetDescription();
	ArgumentFlag GetFlag();
private:
	T* handle;
	void SetToValue(T value);
};

class ArgumentListInterface {
public:
	virtual ~ArgumentListInterface() {}
	virtual ArgumentWrapperInterface& operator[](unsigned index) = 0;
	virtual std::string GetTypeName() = 0;
};

template <typename T>
class ArgumentList : public ArgumentListInterface {
public:
	ArgumentList();
	unsigned int Size() const { return arglist.size(); }
	void Add(ArgumentWrapper<T> arg) { arglist.push_back(arg); }
	ArgumentWrapperInterface& operator[](unsigned index) { return arglist[index]; }
	std::string GetTypeName() { return typestring; }
private:
	std::vector <ArgumentWrapper<T> > arglist;
	std::string typestring;
};

struct MapIndex {
	MapIndex(ArgumentListInterface& arglist, int index) : arglist(arglist), index(index) { }
	ArgumentListInterface& arglist;
	int index;
};

class ArgumentMap {
public:
	template <typename T> void Add(std::string name, ArgumentWrapper<T> variable);
	ArgumentWrapperInterface& Get(std::string name);
	std::string DescriptionBuilder(std::string arg);
private:
	std::map <std::string, MapIndex> argMap;
	
	template <typename T> ArgumentList<T>& ListSelector();
	ArgumentListInterface& GetList(std::string name);
	
	ArgumentList <double> argListDouble;
	ArgumentList <int> argListInt;
	ArgumentList <std::string> argListString;
	ArgumentList <bool> argListBool;
	ArgumentList <unsigned> argListUnsigned;
};

class ArgumentHandler {
public:
	void Parse(int argc, char* argv[]);
	void PrintUsageSyntax();

	template <typename T> void AddRequiredArgument(std::string name, T& variable, std::string description);
	template <typename T, typename U> void AddOptionalArgument(std::string name, T& variable, U defaultValue, std::string description);
	template <typename T, typename U> void AddFlexibleSwitch(std::string name, T& variable, U switchOnValue, U defaultValue, std::string description);
	void AddSwitch(std::string name, bool& variable, std::string description);
private:
	struct BoolString {
		BoolString() : set(false) {}
		BoolString(std::string str) : strval(str), set(true) {}
		BoolString& operator= (std::string str) { strval = str; set = true; return *this; }
		BoolString& operator= (bool flag) { strval = ""; set = flag; return *this; }
		operator std::string() const { return strval; }
		operator bool() const { return set; }
	private:
		std::string strval;
		bool set;
	};
	ArgumentMap expectedArgs;
	std::vector <std::string> argNames;
	std::string programName;
	std::map <std::string, BoolString> readArgs;

	void ParsingReader(int argc, char* argv[]);
};

template <typename T>
void ArgumentHandler::AddRequiredArgument(std::string name, T& variable, std::string description) {
	ArgumentWrapper<T> temp(variable);
	temp.flag = REQUIRED_ARGUMENT;
	temp.description = description;
	expectedArgs.Add(name, temp);
	argNames.push_back(name);
}

template <typename T, typename U>
void ArgumentHandler::AddOptionalArgument(std::string name, T& variable, U defaultValue, std::string description) {
	ArgumentWrapper<T> temp(variable);
	temp.flag = OPTIONAL_ARGUMENT;
	temp.defaultValue = defaultValue;
	temp.description = description;
	expectedArgs.Add(name, temp);
	argNames.push_back(name);
}

template <typename T, typename U>
void ArgumentHandler::AddFlexibleSwitch(std::string name, T& variable, U switchOnValue, U defaultValue, std::string description) {
	ArgumentWrapper<T> temp(variable);
	temp.flag = FLEXIBLE_SWITCH;
	temp.defaultValue = defaultValue;
	temp.switchValue = switchOnValue;
	temp.description = description;
	expectedArgs.Add(name, temp);
	argNames.push_back(name);
}

#endif
