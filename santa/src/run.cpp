#include "run.hpp"

void run_mpi(const char* filename, MPIIO &hIO, int workerNum, int memSize, int lenBuf, unsigned int seed, std::vector<long double > * kerns, double &srcCompCost, double &workerCompCostMax, double &workerCompCostSum)
{
    //Computational Cost
    clock_t begin = clock();

    hIO.init(lenBuf, workerNum);

    // Sourcce init
    if (hIO.isMaster())
    {
        Source source;
        Edge edge;
        EdgeParser parser(filename);
        MID dst1(0);
        MID dst2(0);
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

        // send special edge to mark halfway
        edge.src = numeric_limits<VID>::max(); 
        edge.dst = numeric_limits<VID>::max();
        hIO.bCastEdge(edge, dst1, dst2);

        //cout << " master half" << endl;


        // compute traces of normalize laplacians
        long double termEdges = 0;
        long double termEdgesSq = 0;
        long double termTriangles;
        long double termTwoPaths;
        long double termFourCycles;

        // SECOND PASS
        EdgeParser parserSecondPass(filename);

        while (parserSecondPass.getEdge(edge)) // Stream edges
        {
            if(edge.src != edge.dst) {
                long double currEdgeVal = 1 / ( nodeToDeg[edge.src] * nodeToDeg[edge.dst] );
                termEdges += currEdgeVal;
                termEdgesSq += currEdgeVal*currEdgeVal;
                hIO.bCastEdge(edge, dst1, dst2);
            }
        }

        hIO.sendEndSignal();

        //cout << " master end" << endl;

        // Gather results from curWorkers
        std::vector<long double> workerCnt;
        workerCnt.resize(3,0);
        
        // communication cost for gather
        hIO.recvCnt(source.getMaxVId(), workerCnt);
        // final node count
        n++;

        hIO.recvTime(workerCompCostMax, workerCompCostSum);
        
        for (int i = 0; i < 3; i++){
            workerCnt[i] = workerCnt[i] / workerNum;
        }

        termTwoPaths = workerCnt[0];
        termTriangles = workerCnt[1];
        termFourCycles = workerCnt[2];

        // compute traces
        long double traceL1 = n;
        long double traceL2 = n + 2*termEdges;
        long double traceL3 = n + 6*termEdges + 6*termTriangles;
        long double traceL4 = n + 12*termEdges + 2*termEdgesSq + 
                                        24*termTriangles + 4*termTwoPaths + 8*termFourCycles;

        // create descriptor
        for ( int i = 0; i < DESC_SIZE ; i++ ){
            long double t = arrayOfTs[i];
            long double term0 = n;
            long double term1 = t * traceL1;
            long double term2 = t * t * traceL2 / 2;
            long double term3 = t * t * t * traceL3 / 6;
            long double term4 = t * t * t * t * traceL4 / 24;

            long double heatSum = term0 - term1 + term2 - term3 + term4;
            long double waveSum = term0 - term2 + term4;

            kerns[0][i] = heatSum; // heat
            kerns[1][i] = heatSum / n; // heat, empty normalization 
            kerns[2][i] = heatSum / ( 1 + ( n - 1 )*exp( -1 * t ) ); // heat, complete normalization
            kerns[3][i] = waveSum; // wave
            kerns[4][i] = waveSum / n; // wave, empty normalization
            kerns[5][i] = waveSum / ( 1 + ( n - 1 )*cos(t) ); // wave, complete normalization
        }

        srcCompCost = (double(clock() - begin)) / CLOCKS_PER_SEC; // source cpu time
        
        //std::cout << elapsedTime1 << "\t" << elapsedTime  << endl;
        //std::cout << hIO.getCommCostDistribute() << endl;
        //std::cout << srcCompCost  << endl;
        //std::cout << "master ends..." << endl;

        // report results
        
        //cout << " master total end " << endl;

        return;
    }
    else // Worker part
    {

        //std::cout << "worker begins..." << endl;

        Worker  worker(memSize, seed + hIO.getWorkerId());
        Edge edge;

        while(hIO.recvEdge(edge))
        {
            if ( edge.src == numeric_limits<VID>::max() && edge.dst == numeric_limits<VID>::max() )
                break;
            worker.processEdge1(edge);
        }
        //std::cout << "worker 222 begins... " << endl;

        // second pass
        while(hIO.recvEdge(edge))
        {
            worker.processEdge2(edge);
        }

        //std::cout << "Worker: "  << double(clock() - begin) / CLOCKS_PER_SEC << "\t" << hIO.getIOCPUTime() / CLOCKS_PER_SEC << "\t" <<  workerCompCost << endl;
        //std::cout << "worker ends" << endl;

        // send counts to master
        
        hIO.sendCnt(worker.getTerms());

        double workerCompCost = (double(clock() - begin) - hIO.getIOCPUTime()) / CLOCKS_PER_SEC; // source cpu time

        hIO.sendTime(workerCompCost);


        //cout << " worker total end " << endl;
        
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

            //long double kerns[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            std::vector<long double > kerns[6];
            for ( int k = 0; k < 6 ; k++ )
                kerns[k].resize(DESC_SIZE, 0);
            
            run_mpi(input, hIO, workerNum, memSize, bufLen, seed + repeat * workerNum * i, kerns, srcCompCost, workerCompCostMax, workerCompCostSum);
            
            gettimeofday(&endTV, NULL);

            timersub(&endTV, &startTV, &diff);

            //double elapsedTime = diff.tv_sec * 1000 + diff.tv_usec / 1000 ;
            
            for ( int k = 0; k < 6 ; k++ )
                print_cnt(outPath, kerns, i, k);

            //cout << " print time start " << endl;

            std::ostringstream timeFileName;
            timeFileName << outPath << "/time" << i <<".txt";
            std::fstream tfp;
            tfp.open(timeFileName.str(), std::fstream::out | std::fstream::trunc);
            tfp << srcCompCost << endl;
            tfp.close();
            
            //cout << " print time end " << endl;

        } else {

            double srcCompCost = 0;
            double workerCompCostMax = 0;
            double workerCompCostSum = 0;
            std::vector<long double > kerns[6];
            run_mpi(input, hIO, workerNum, memSize, bufLen, seed + repeat * workerNum * i, kerns, srcCompCost, workerCompCostMax, workerCompCostSum);
            
        }
    }
}


void print_cnt(const char* outPath, const std::vector<long double> * kerns, int id, int type)
{
	// Print global count
	std::ostringstream tfileName;
	tfileName << outPath;

    switch( type ){
        case 0:
            tfileName << "/heat";
            break;
        case 1:
            tfileName << "/heatEmpty";
            break;
        case 2:
            tfileName << "/heatComplete";
            break;
        case 3:
            tfileName << "/wave";
            break;
        case 4:
            tfileName << "/waveEmpty";
            break;
        default:
            tfileName << "/waveComplete";
            break;
    }

    tfileName << id << ".txt";
	std::fstream	gfp;
	gfp.open(tfileName.str(), std::fstream::out | std::fstream::trunc);
    for (int i = 0; i < DESC_SIZE; i++){
        gfp << std::setprecision(std::numeric_limits<double>::max_digits10) << kerns[type][i] << endl;
    }
	gfp.close();

        //cout << " print total end " << endl;
}