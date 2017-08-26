#include "libraries/operaArgumentHandler.h"

operaArgumentHandler::operaArgumentHandler() {
	AddSwitch("verbose", verbose, "Output informational messages");
	AddSwitch("debug", debug, "Output debug messages");
	AddSwitch("trace", trace, "Output trace messages");
	AddSwitch("plot", plot, "Produce plots");
}

void operaArgumentHandler::AddPlotFileArguments(std::string& plotfilename, std::string& datafilename, std::string& scriptfilename, bool& interactive) {
	AddOptionalArgument("plotfilename", plotfilename, "", "Output plot eps file name");
	AddOptionalArgument("datafilename", datafilename, "", "Output data file name");
	AddOptionalArgument("scriptfilename", scriptfilename, "", "Output gnuplot script file name");
	AddSwitch("interactive", interactive, "For interactive plots");
}

void operaArgumentHandler::AddOrderLimitArguments(int& ordernumber, int& minorder, int& maxorder, const int default_value) {
	AddOptionalArgument("ordernumber", ordernumber, default_value, "Specific single order to use");
	AddOptionalArgument("minorder", minorder, default_value, "Minimum order to use");
	AddOptionalArgument("maxorder", maxorder, default_value, "Maximum order to use");
}

bool operaArgumentHandler::verbose = false;
bool operaArgumentHandler::debug = false;
bool operaArgumentHandler::trace = false;
bool operaArgumentHandler::plot = false;
