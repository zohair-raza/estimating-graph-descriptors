#include "run.hpp"

void run_mpi(const char* filename, MPIIO &hIO, int workerNum, int memSize, int lenBuf, unsigned int seed, std::vector<long double > & globalCnt, double &srcCompCost, double &workerCompCostMax, double &workerCompCostSum)
{

    //Computational Cost
    clock_t begin = clock();

    hIO.init(lenBuf, workerNum);

    // Sourcce init
    if (hIO.isMaster())
    {
        Source source;
        Edge edge;
        MID dst1(0);
        MID dst2(0);
        EdgeParser parser(filename);
        long double m = 0;
        long double n = 0;
        
        std::vector<long double> nodeToDeg;
        nodeToDeg.resize(1000000, 0);
        
        
        //cout << " master start" << endl;
        
        while (parser.getEdge(edge)) // Stream edges
        {
            if(edge.src != edge.dst) {
                m++;
                
                if (edge.src > n) {
                    n = edge.src;
                }
                if (edge.dst > n) {
                    n = edge.dst;
                }
                if (n > nodeToDeg.size()) {
                    nodeToDeg.resize(2*n,0);
                }
                
                //cout << " master sending edge "<< edge.src << "-" << edge.dst << endl;
                
                
                globalCnt[4] += nodeToDeg[edge.src] + nodeToDeg[edge.dst];
                globalCnt[11] += (nodeToDeg[edge.src]*(nodeToDeg[edge.src]-1))/2 + (nodeToDeg[edge.dst]*(nodeToDeg[edge.dst]-1))/2;
                
                //cout << " master sending edge 222"<< edge.src << "-" << edge.dst << endl;
                
                nodeToDeg[edge.src] += 1;
                nodeToDeg[edge.dst] += 1;
                
                
                source.processEdge(edge, dst1, dst2);
                hIO.bCastEdge(edge, dst1, dst2);
            }
        }

        //cout << " master done 0" << endl;
        
        hIO.sendEndSignal();

        //std::cout << "Master: " << double(clock() - begin) / CLOCKS_PER_SEC << "\t" << hIO.getIOCPUTime() / CLOCKS_PER_SEC << "\t" <<  srcCompCost << endl;

        // Gather results from curWorkers
        
        std::vector<long double> workerCnt;
        workerCnt.resize(6,0);
        
        // communication cost for gather
        hIO.recvCnt(source.getMaxVId(), workerCnt);
        n++;
        
        //std::cout << source.getMaxVId() << "\t" << globalCnt << "\t" << oLocalCnt.size();

        hIO.recvTime(workerCompCostMax, workerCompCostSum);
        
        for (int i = 0; i < 6; i++){
            workerCnt[i] = workerCnt[i] / workerNum;
        }
        
        
        long double nc2 = (n*(n-1))/2;
        /*long double nc3 = nc2*(n-2)/3;
        long double nc4 = nc3*(n-3)/4;
        
        // motif counts
        
        globalCnt[0] = n*(n-1)/2;
        globalCnt[1] = m;
        
        globalCnt[2] = n*(n-1)*(n-2)/6;
        globalCnt[3] = m*(n-2);
        globalCnt[5] = workerCnt[0];
        
        globalCnt[6] = globalCnt[2]*(n-3)/4;
        globalCnt[7] = m/2*(n-2)*(n-3);
        globalCnt[8] = globalCnt[4]*(n-3);
        globalCnt[9] = (m*(m-1)/2) - globalCnt[4];
        globalCnt[10] = globalCnt[5]*(n-3);
        
        globalCnt[12] = workerCnt[1];
        globalCnt[13] = workerCnt[2];
        globalCnt[14] = workerCnt[3];
        globalCnt[15] = workerCnt[4];
        globalCnt[16] = workerCnt[5];
        */
        
        // pre-normalized motif counts
        globalCnt[0] = 1;
        globalCnt[1] = m/(nc2);
        
        globalCnt[2] = 1;
        globalCnt[3] = globalCnt[1]*3;
        globalCnt[4] = (((globalCnt[4]/n)/(n-1))*6)/(n-2);
        globalCnt[5] = (((workerCnt[0]/n)/(n-1))*6)/(n-2);
        
        globalCnt[6] = 1;
        globalCnt[7] = globalCnt[1]*6;
        globalCnt[8] = globalCnt[4]*4;
        long double temp = (m-1)*6/((n-2)*(n-3));
        globalCnt[9] = globalCnt[1]*temp - (globalCnt[4]*4)/(n-3);
        globalCnt[10] = globalCnt[5]*4;
        globalCnt[11] = ((((globalCnt[11]/n)/(n-1))*24)/(n-2)/(n-3));
        
        globalCnt[12] = ((((workerCnt[1]/n)/(n-1))*24)/(n-2)/(n-3));
        globalCnt[13] = ((((workerCnt[2]/n)/(n-1))*24)/(n-2)/(n-3));
        globalCnt[14] = ((((workerCnt[3]/n)/(n-1))*24)/(n-2)/(n-3));
        globalCnt[15] = ((((workerCnt[4]/n)/(n-1))*24)/(n-2)/(n-3));
        globalCnt[16] = ((((workerCnt[5]/n)/(n-1))*24)/(n-2)/(n-3));
        
        
        // graphlet counts
        
        globalCnt[0] = globalCnt[0] - globalCnt[1];
        
        globalCnt[2] = globalCnt[2] - globalCnt[3] + globalCnt[4]- globalCnt[5];
        globalCnt[3] = globalCnt[3] - 2*globalCnt[4] + 3*globalCnt[5];
        globalCnt[4] = globalCnt[4] - 3*globalCnt[5];
        
        globalCnt[6] = globalCnt[6] - globalCnt[7] + globalCnt[8] + globalCnt[9] - globalCnt[10] - globalCnt[11] - globalCnt[12] + globalCnt[13] + globalCnt[14] - globalCnt[15] + globalCnt[16];
        globalCnt[7] = globalCnt[7] -2*globalCnt[8] - 2*globalCnt[9] + 3*globalCnt[10] + 3*globalCnt[11] + 3*globalCnt[12] - 4*globalCnt[13] - 4*globalCnt[14] + 5*globalCnt[15] - 6*globalCnt[16];
        globalCnt[8] = globalCnt[8] -3*globalCnt[10]-3*globalCnt[11]-2*globalCnt[12]+ 5*globalCnt[13]+ 4*globalCnt[14]-8*globalCnt[15]+ 12*globalCnt[16];
        globalCnt[9] = globalCnt[9] -globalCnt[12]+ globalCnt[13]+ 2*globalCnt[14]- 2*globalCnt[15]+ 3*globalCnt[16];
        globalCnt[10] = globalCnt[10] - globalCnt[13]  + 2*globalCnt[15] -4*globalCnt[16];
        globalCnt[11] = globalCnt[11] - globalCnt[13]  + 2*globalCnt[15] -4*globalCnt[16];
        globalCnt[12] = globalCnt[12] - 2*globalCnt[13] -4*globalCnt[14] + 6*globalCnt[15] - 12*globalCnt[16];
        globalCnt[13] = globalCnt[13] - 4* globalCnt[15] + 12*globalCnt[16];
        globalCnt[14] = globalCnt[14] - globalCnt[15] + 3*globalCnt[16];
        globalCnt[15] = globalCnt[15] - 6*globalCnt[16];
        
        
        // normalize
        /*
        for (int i = 0; i < 2; i++){
            globalCnt[i] /= nc2;
        }
        for (int i = 2; i < 6; i++){
            globalCnt[i] /= nc3;
        }
        for (int i = 6; i < 17; i++){
            globalCnt[i] /= nc4;
        }
        */
        
        srcCompCost = (double(clock() - begin)) / CLOCKS_PER_SEC; // source cpu time
        
        //std::cout << elapsedTime1 << "\t" << elapsedTime  << endl;
        //std::cout << hIO.getCommCostDistribute() << endl;
        //std::cout << srcCompCost  << endl;
        //std::cout << "master ends..." << endl;

        // report results
        
        return;
    }
    else // Worker part
    {

        //std::cout << "worker begins..." << endl;

        Worker  worker(memSize, seed + hIO.getWorkerId());
        Edge edge;
        while(hIO.recvEdge(edge))
        {
            worker.processEdge(edge);
        }

        //std::cout << "Worker: "  << double(clock() - begin) / CLOCKS_PER_SEC << "\t" << hIO.getIOCPUTime() / CLOCKS_PER_SEC << "\t" <<  workerCompCost << endl;
        //std::cout << "worker ends" << endl;

        // send counts to master
        
        
        hIO.sendCnt(worker.getGlobalCnt());

        double workerCompCost = (double(clock() - begin) - hIO.getIOCPUTime()) / CLOCKS_PER_SEC; // source cpu time

        hIO.sendTime(workerCompCost);
        
        return;
    }

}

void run_exp (const char* input, const char* outPath, MPIIO &hIO, int workerNum, int memSize, int repeat, int bufLen)
{

    int seed = 0;

    struct timeval diff, startTV, endTV;

	if (hIO.isMaster())
	{
		struct stat sb;
		if (stat(outPath, &sb) == 0)
		{
			if (S_ISDIR(sb.st_mode)) //TODO. directory is exists
				;
			else if (S_ISREG(sb.st_mode)) //TODO. No directory but a regular file with same name
				;
			else // TODO. handle undefined cases.
				;
		} 
		else 
		{
			mkdir(outPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}
    }

    for(int i =0 ; i < repeat; i++) {

        if (hIO.isMaster()) {

            gettimeofday(&startTV, NULL);

            double srcCompCost = 0;
            double workerCompCostMax = 0;
            double workerCompCostSum = 0;

            //long double globalCnt[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            std::vector<long double > globalCnt;
            globalCnt.resize(17,0);
            
            
            run_mpi(input, hIO, workerNum, memSize, bufLen, seed + repeat * workerNum * i, globalCnt, srcCompCost, workerCompCostMax, workerCompCostSum);
            
            gettimeofday(&endTV, NULL);

            timersub(&endTV, &startTV, &diff);

            //double elapsedTime = diff.tv_sec * 1000 + diff.tv_usec / 1000 ;
            
            print_cnt(outPath, globalCnt, i);
            std::ostringstream timeFileName;
            timeFileName << outPath << "/time" << i <<".txt";
            std::fstream tfp;
            tfp.open(timeFileName.str(), std::fstream::out | std::fstream::trunc);
            tfp << srcCompCost << endl;
            tfp.close();
            
        } else {

            double srcCompCost = 0;
            double workerCompCostMax = 0;
            double workerCompCostSum = 0;
            std::vector<long double > globalCnt;
            run_mpi(input, hIO, workerNum, memSize, bufLen, seed + repeat * workerNum * i, globalCnt, srcCompCost, workerCompCostMax, workerCompCostSum);
            
        }
    }
}

void print_cnt(const char* outPath, const std::vector<long double> &globalCnt, int id)
{

	// Print global count
	std::ostringstream gCntFileName;
	gCntFileName << outPath << "/global" << id << ".txt";
	std::fstream	gfp;
	gfp.open(gCntFileName.str(), std::fstream::out | std::fstream::trunc);
    for (int i = 0; i < 17; i++){
        gfp << std::setprecision(std::numeric_limits<double>::max_digits10) << globalCnt[i] << endl;
    }
	gfp.close();
}
