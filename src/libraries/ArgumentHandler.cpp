#include "libraries/ArgumentHandler.h"
#include <iostream>
#include <stdexcept>
#include <sstream>

template <typename T>
void ArgumentWrapper<T>::SetToDefaultValue() {
	if (flag == REQUIRED_ARGUMENT) throw std::runtime_error("missing required argument");
	SetValue(defaultValue);
}

template <typename T>
void ArgumentWrapper<T>::SetToSwitchValue() {
	if (flag == REQUIRED_ARGUMENT || flag == OPTIONAL_ARGUMENT) throw std::runtime_error("argument expects a value");
	SetValue(switchValue);
}

template <typename T>
void ArgumentWrapper<T>::SetToValue(T value) {
	if (flag == SWITCH) throw std::runtime_error("argument does not expect a value");
	SetValue(value);
}

template <typename T>
void ArgumentWrapper<T>::SetToValueFromString(std::string value) {
	T temp;
	std::istringstream ss(value);
	ss >> temp; //future work: use more stringent type-checking and error throwing?
	if (ss.fail()) throw std::runtime_error("value could not be read into argument");
	SetToValue(temp);
}

template <>
void ArgumentWrapper<std::string>::SetToValueFromString(std::string value) {
	SetToValue(value);
}

template <typename T>
std::string ArgumentWrapper<T>::GetDescription() {
	return description;
}

template <typename T>
ArgumentFlag ArgumentWrapper<T>::GetFlag() {
	return flag;
}

template class ArgumentWrapper<double>;
template class ArgumentWrapper<int>;
template class ArgumentWrapper<std::string>;
template class ArgumentWrapper<bool>;
template class ArgumentWrapper<unsigned>;

template <> ArgumentList<double>::ArgumentList() : typestring("double") { }
template <> ArgumentList<int>::ArgumentList() : typestring("int") { }
template <> ArgumentList<std::string>::ArgumentList() : typestring("string") { }
template <> ArgumentList<bool>::ArgumentList() : typestring("bool") { }
template <> ArgumentList<unsigned>::ArgumentList() : typestring("unsigned") { }

template <> ArgumentList<double>& ArgumentMap::ListSelector<double>() { return argListDouble; }
template <> ArgumentList<int>& ArgumentMap::ListSelector<int>() { return argListInt; }
template <> ArgumentList<std::string>& ArgumentMap::ListSelector<std::string>() { return argListString; }
template <> ArgumentList<bool>& ArgumentMap::ListSelector<bool>() { return argListBool; }
template <> ArgumentList<unsigned>& ArgumentMap::ListSelector<unsigned>() { return argListUnsigned; }

template void ArgumentMap::Add<double>(std::string name, ArgumentWrapper<double> variable);
template void ArgumentMap::Add<int>(std::string name, ArgumentWrapper<int> variable);
template void ArgumentMap::Add<std::string>(std::string name, ArgumentWrapper<std::string> variable);
template void ArgumentMap::Add<bool>(std::string name, ArgumentWrapper<bool> variable);
template void ArgumentMap::Add<unsigned>(std::string name, ArgumentWrapper<unsigned> variable);

std::string ArgumentMap::DescriptionBuilder(std::string arg) {
	ArgumentWrapperInterface& argref = Get(arg);
	ArgumentFlag flag = argref.GetFlag();
	std::string temp;
	if (flag != REQUIRED_ARGUMENT) temp += "[";
	temp += "--";
	temp += arg;
	if (flag != SWITCH) {
		if (flag == FLEXIBLE_SWITCH) temp += "[";
		temp += "=";
		temp += GetList(arg).GetTypeName();
		if (flag == FLEXIBLE_SWITCH) temp += "]";
	}
	if (flag != REQUIRED_ARGUMENT) temp += "]";
	temp += "\t";
	temp += argref.GetDescription();
	return temp;
}

template <typename T>
void ArgumentMap::Add(std::string name, ArgumentWrapper<T> variable) {
	ArgumentList<T>& argListTyped = dynamic_cast<ArgumentList<T>&>(ListSelector<T>());
	argListTyped.Add(variable);
	if (argMap.count(name) != 0) throw std::logic_error("multiple arguments mapped to the name " + name);
	MapIndex ind(argListTyped, argListTyped.Size() - 1);
	argMap.insert(std::make_pair(name, ind));
}

ArgumentWrapperInterface& ArgumentMap::Get(std::string name) {
	std::map <std::string, MapIndex>::const_iterator temp = argMap.find(name);
	if (temp == argMap.end()) throw std::logic_error("no argument in argument map named " + name);
	return temp->second.arglist[temp->second.index];
}

ArgumentListInterface& ArgumentMap::GetList(std::string name) {
	std::map <std::string, MapIndex>::const_iterator temp = argMap.find(name);
	if (temp == argMap.end()) throw std::logic_error("no argument in argument map named " + name);
	return temp->second.arglist;
}

void ArgumentHandler::AddSwitch(std::string name, bool& variable, std::string description) {
	ArgumentWrapper<bool> temp(variable);
	temp.flag = SWITCH;
	temp.defaultValue = false;
	temp.switchValue = true;
	temp.description = description;
	expectedArgs.Add(name, temp);
	argNames.push_back(name);
}

void ArgumentHandler::ParsingReader(int argc, char* argv[]) {
	programName = argv[0];
	for (int i = 1; i < argc; i++) {
		std::string temp = argv[i];
		if (temp.length() < 3 || temp[0] != '-' || temp[1] != '-') throw std::runtime_error("expected -- followed by argument name, got " + temp);
		std::size_t splitpos = temp.find('=');
		std::string readName = temp.substr(2, splitpos - 2);
		BoolString readVal = temp.substr(splitpos + 1, std::string::npos);
		if (splitpos == std::string::npos) readVal = false;
		if (readArgs.count(readName) > 0) throw std::runtime_error("same argument name passed multiple times: " + readName);
		readArgs[readName] = readVal;
	}
}

void ArgumentHandler::Parse(int argc, char* argv[]) {
	class HelpMessage : std::exception { };
	try {
		ParsingReader(argc, argv);
		if (readArgs.count("help") > 0) {
			throw HelpMessage();
		}
		for (unsigned int i = 0; i < argNames.size(); i++) {
			try {
				ArgumentWrapperInterface& argref = expectedArgs.Get(argNames[i]);
				if (readArgs.count(argNames[i]) == 0) argref.SetToDefaultValue();
				else if (!readArgs[argNames[i]]) argref.SetToSwitchValue();
				else argref.SetToValueFromString(readArgs[argNames[i]]);
			}
			catch (std::runtime_error& error) {
				throw std::runtime_error(error.what() + (": " + argNames[i]));
			}
			readArgs.erase(argNames[i]);
		}
		if (readArgs.size() > 0) throw std::runtime_error("unknown argument read in: " + readArgs.begin()->first);
	}
	catch (HelpMessage) {
		PrintUsageSyntax();
		throw std::exception();
	}
	catch (std::exception& error) {
		std::cout << "Error: " << error.what() << std::endl;
		std::cout << "See --help for available options." << std::endl;
		throw error;
	}
}

void ArgumentHandler::PrintUsageSyntax() {
	std::cout << "Usage options" << (programName.empty() ? "" : " for ") << programName << ":" << std::endl;
	for (unsigned int i = 0; i < argNames.size(); i++) {
		std::cout << expectedArgs.DescriptionBuilder(argNames[i]) << std::endl;
	}
}
