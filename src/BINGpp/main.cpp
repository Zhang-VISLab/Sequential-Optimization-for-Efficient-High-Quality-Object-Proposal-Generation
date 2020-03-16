//This source code is free for academic usage under LICENSE.

#include "stdafx.h"
#include "BINGPP.h"
#include "ValStructVec.h"
#include "CmShow.h"

void RunBINGPP(CStr &resName, double base, int W, int NSS, int numPerSz, CStr &dataPath);

int main(int argc, char* argv[])
{
	//DataSetVOC::importImageNetBenchMark();
	//DataSetVOC::cvt2OpenCVYml("C:/WkDir/DetectionProposals/VOC2007/Annotations/");
	if(argc < 2){
		std::cerr << "Please pass the data path to as first argument" << std::endl;
		return 1;
	}
	CStr dataPath = argv[1];
	RunBINGPP("WinRecall.m", 2, 8, 2, 130, dataPath);
	return 0;
}

void RunBINGPP(CStr &resName, double base, int W, int NSS, int numPerSz, CStr &dataPath)
{
	srand(131);//srand((unsigned int)time(NULL));
	//omp_set_num_threads(16);
	DataSetVOC voc2007(dataPath); 
	voc2007.loadAnnotations();
	//voc2007.loadDataGenericOverCls();

	cout << "Dataset:'" << _S(voc2007.wkDir) << "' with " << voc2007.trainNum << " training and " << voc2007.testNum << " testing" << endl;
	cout << _S(resName) << " Base = " << base << ", W = " << W << ", NSS = " << NSS << ", perSz = " << numPerSz << endl;

	BINGPP BINGpp(voc2007, base, W, NSS);
    
	vector<vector<Vec4i>> boxes;
	//objNess.getObjBndBoxesForTests(boxes, 250);
         
	BINGpp.getObjBndBoxesForTestFast(boxes, numPerSz);
	//objNess.getObjBndBoxesForTrainsEva(boxes, numPerSz);
	//objNess.getRandomBoxes(boxes);
	//objNess.evaluatePerClassRecall(boxes, resName, 2000);
	//objNess.illuTestReults(boxes);
}

