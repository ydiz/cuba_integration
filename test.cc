// #include <qlat/qlat.h>

#include "integration.h"
#include "integration_init.h"

int main(int argc, char **argv) {

	Int_para para;
	init(argc, argv, para);

	if(para.term == "all_with_p1") run_p1(para);	
	else run(para);
	return 0;
}
