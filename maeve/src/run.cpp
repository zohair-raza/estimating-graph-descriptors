#include "run.hpp"
#include <math.h>






void run_mpi(const char* filename, MPIIO &hIO, int workerNum, int memSize, int lenBuf, unsigned int seed, long double * moments, double &srcCompCost, double &workerCompCostMax, double &workerCompCostSum)
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
                
                nodeToDeg[edge.src] += 1;
                nodeToDeg[edge.dst] += 1;
                
                source.processEdge(edge, dst1, dst2);
                hIO.bCastEdge(edge, dst1, dst2);
            }
        }

        hIO.sendEndSignal();

        //std::cout << "Master: " << double(clock() - begin) / CLOCKS_PER_SEC << "\t" << hIO.getIOCPUTime() / CLOCKS_PER_SEC << "\t" <<  srcCompCost << endl;
        
        
        std::vector<long double> oLocalCnt[2];
        
        // Gather results from curWorkers
        hIO.recvCnt(source.getMaxVId(), oLocalCnt[0]);
        hIO.recvCnt(source.getMaxVId(), oLocalCnt[1]);

        // communication cost for gather
        
        n++;
        
        hIO.recvTime(workerCompCostMax, workerCompCostSum);
        
        for (int i = 0; i < 2; i++){
            for(auto it = oLocalCnt[i].begin(); it != oLocalCnt[i].end(); ++it)
            {
                *it  = *it / workerNum;
            }
        }
        
        //std::vector<long double> feature1;
        //feature1.resize(n, 0);
        
        //std::vector<long double> moments;
        //moments->resize(20,0);
        
        long double avg[5] = {0.0,0.0,0.0,0.0,0.0};
        long double var[5] = {0.0,0.0,0.0,0.0,0.0};
        long double skew[5] = {0.0,0.0,0.0,0.0,0.0};
        long double kurt[5] = {0.0,0.0,0.0,0.0,0.0};
        
        
        avg[0] = 2*m/n;
        
        for (unsigned int i = 0; i < n; i++){
            
            if (nodeToDeg[i] == 0){
                continue;
            }
            
            long double f1 = (nodeToDeg[i] == 1) ? 0 : 2*oLocalCnt[0][i]/(nodeToDeg[i]*(nodeToDeg[i]-1));
            
            avg[1] += (f1)/n;
            avg[2] += (1 + oLocalCnt[1][i]/nodeToDeg[i] )/n;
            avg[3] += (nodeToDeg[i] + oLocalCnt[0][i])/n;
            avg[4] += (oLocalCnt[1][i] - 2*oLocalCnt[0][i])/n;
            
        }
        
        for (unsigned int i = 0; i < n; i++){
            
            
            long double f[5] = {nodeToDeg[i], (nodeToDeg[i] == 1) ? 0 : 2*oLocalCnt[0][i]/(nodeToDeg[i]*(nodeToDeg[i]-1)), 1 + oLocalCnt[1][i]/nodeToDeg[i], nodeToDeg[i] + oLocalCnt[0][i], oLocalCnt[1][i] - 2*oLocalCnt[0][i]};
            
            if (nodeToDeg[i] == 0){
                for (int j = 1; j < 5; j++){
                    f[j] = 0;
                }
            }
            
            for (int j = 0; j < 5; j++){
                var[j] += pow(f[j]- avg[j], 2)/n;
                skew[j] += pow(f[j]- avg[j], 3)/n;
                kurt[j] += pow(f[j]- avg[j], 4)/n;
            }
        }
        
        long double stdev[5] = {0.0,0.0,0.0,0.0,0.0};
        for (int i = 0; i < 5; i++){
            stdev[i] = sqrt(var[i]);
            if (var[i] != 0){
                skew[i] /= pow(stdev[i], 3);
                kurt[i] /= pow(var[i],2);
            }
        }
        
        
        for (int i = 0; i < 5; i++){
            moments[4*i] = avg[i];
            moments[4*i+1] = stdev[i];
            moments[4*i+2] = skew[i];
            moments[4*i+3] = kurt[i];
        }
        
        //srcCompCost = (double(clock() - begin) - hIO.getIOCPUTime()) / CLOCKS_PER_SEC; // source cpu time
        
        srcCompCost = (double(clock() - begin))/ CLOCKS_PER_SEC;
        
        //cout << " time : " << srcCompCost << endl;
        
        
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
        hIO.sendCnt(worker.getLocalCnt(0));
        hIO.sendCnt(worker.getLocalCnt(1));

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

            long double moments[20];

            double srcCompCost = 0;
            double workerCompCostMax = 0;
            double workerCompCostSum = 0;

            
            run_mpi(input, hIO, workerNum, memSize, bufLen, seed + repeat * workerNum * i, moments, srcCompCost, workerCompCostMax, workerCompCostSum);

            
            gettimeofday(&endTV, NULL);

            timersub(&endTV, &startTV, &diff);

            //double elapsedTime = diff.tv_sec * 1000 + diff.tv_usec / 1000 ;
            
            print_cnt(outPath, moments, i);

            // print time
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
            long double moments[20];
            run_mpi(input, hIO, workerNum, memSize, bufLen, seed + repeat * workerNum * i, moments, srcCompCost, workerCompCostMax, workerCompCostSum);
        }
    }
}

void print_cnt(const char* outPath, const long double * moments, int id)
{

	// Print global count

	std::ostringstream lCntFileName;
	lCntFileName << outPath << "/local" << id <<".txt";
	std::fstream	lfp;
	lfp.open(lCntFileName.str(), std::fstream::out | std::fstream::trunc);

	for (unsigned int nid = 0; nid < 20; nid++)
	{
		lfp << moments[nid] << endl;
	}
	lfp.close();
}
