#include "integration.h"
#include "integration_init.h"

int main(int argc, char **argv) {

	Int_para para;
	init(argc, argv, para);

	if(para.term == "with_p1") run_p1(para);	
	else if(para.term == "with_p3")  run_p3(para);
  else assert(0);
	return 0;
}
