#include "FoldingApplication.h"


int main(int argc, char* argv[])
{
	std::string windowNameArgv;
	for (int i = 0; i < argc; ++i) {
		std::string str(argv[i]);
		windowNameArgv += " ";
		windowNameArgv += str;
	}
	char* a[] = { windowNameArgv.data() };
	FoldingApplication app(1, a);
	app.launch(argc, argv);
	return 0;

}
