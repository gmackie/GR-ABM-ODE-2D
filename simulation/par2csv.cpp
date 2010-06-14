#include "gr.h"
#include "params.h"

int main(int argc, char** argv)
{
	if (argc != 3)
		return 1;

	bool header = (atoi(argv[1]) != 0);

	if (!Params::getInstance()->fromXml(argv[2]))
		return 1;

	Params::getInstance()->printCSV(header);

	return 0;
}
