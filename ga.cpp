/***************************************************************************
                          ga.cpp  -  MPI with OpenMP algorithm
                          implementation of the original timetabling GA
                             -------------------
    copyright            : (C) 2015 by Neli Lepoeva
    email                : neli.lepoeva@gmail.com
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation - version 3 of the License               *
 *                                                                         *
 ***************************************************************************/
#include "Control.h"
#include "Problem.h"
#include "Solution.h"
#include "Timer.h"

#include "json/json.h"

#include <fstream>

#include <list>
#include <vector>
#include <string.h>
#include <cstring>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>
#include <netdb.h>
#include <errno.h>
#include <arpa/inet.h>
#include <pthread.h>

#include <mpi.h>
#include <omp.h>

using namespace Json;

#define MALLOC_ERROR    -2
#define SEND_RECV_TAG   2
#define SEND_RECV_TAG_2 4

Timer timer;
Random* rnd;
double elapsed_time;

int id; // pocessor ID
int p; //number of processors for the MPI run
int rc; // MPI return code

int t; // number of threads for the OpenMP run
int tid; // thread ID

int bestScv; // best evaluation (with only soft constarint violations)
int bestEvaluation; // best evaluation (there are hard constraint violations)

ostream *os = &cout;
istream *is = &cin;

int problemType;
const int popSize = 10;
Solution* pop[popSize];
int maxSteps; // max number of steps in the local search

Problem *problem;
int sizeOfProblem;
void *bufferProblem;
int sizeOfSolution;
void *bufferSolution;

int seed; // seed value for the random number generator
int n_of_events;
int n_of_rooms;
int n_of_features;
int n_of_students;

/*
 * Returns the problem size in bytes
 *
 * @param The Problem instance
 */
int getProblemSizeinB(Problem *problem) {
    int sizeOf_studentNumbern = problem->n_of_events;
    int sizeOf_roomSize = problem->n_of_rooms;
    int sizeOf_student_events = problem->n_of_students*problem->n_of_events;
    int sizeOf_eventCorrelations = problem->n_of_events*problem->n_of_events;
    int sizeOf_room_features = problem->n_of_rooms*problem->n_of_features;
    int sizeOf_event_features = problem->n_of_events*problem->n_of_features;
    int sizeOf_possibleRooms = problem->n_of_events*problem->n_of_rooms;
    return 4 * sizeof(int) + sizeOf_studentNumbern * sizeof(int) + sizeOf_roomSize * sizeof(int) +
        sizeOf_student_events * sizeof(int) + sizeOf_eventCorrelations * sizeof(int) +
        sizeOf_room_features * sizeof(int) + sizeOf_event_features * sizeof(int) +
        sizeOf_possibleRooms * sizeof(int);
}

/*
 * Returns the solution size in bytes
 * (calculated from the data in the problem)
 *
 * @param The Problem instance
 */
int getSolutionSizeInB(Problem *problem) {
    int sizeOf_solutionsVector = 2 * problem->n_of_events;
    return 3 * sizeof(int) + sizeOf_solutionsVector * sizeof(int) + sizeof(bool);
}

/*
 * Tournament selection with tornament size 2
 */
Solution* selection(Solution** pop, int popSize ){

    // tournament selection with tornament size 2
    int first, second;
    first = (int)(rnd->next()*popSize);
    second = (int)(rnd->next()*popSize);

    if(pop[first]->penalty < pop[second]->penalty)
        return(pop[first]);
    else
        return(pop[second]);
}

/*
 * Tournament selection with tornament size 5
 */
Solution* selection5(Solution** pop, int popSize ) {

    // tournament selection with tornament size 5
    int tournament[5];
    int best;

    tournament[0] = (int)(rnd->next()*popSize);

    best =  tournament[0];
    for(int i= 1; i<5; i++){
        tournament[i] = (int)(rnd->next()*popSize);
        if(pop[tournament[i]]->penalty < pop[best]->penalty)
            best = tournament[i];
    }

    return(pop[best]);
}

/*
 * Compares two solutions by their penalties
 */
bool compareSolution(Solution * sol1, Solution * sol2)
{
    return sol1->penalty < sol2->penalty;
}

void resetTime() {
    timer.resetTime();
}

double getTime() {
    return timer.elapsedTime( Timer::REAL );
}

void beginTry() {
    resetTime();
    bestScv = INT_MAX;
    bestEvaluation = INT_MAX;
}

void endTry(Solution *bestSolution, Value solution) {//, StreamWriterBuilder swb, Value solution) {//
    StreamWriterBuilder swb;
    swb.settings_["indentation"] = "";
    solution["solution"]["threadID"] = tid;
    if (bestSolution->feasible) {
        solution["solution"]["totalTime"] = timer.elapsedTime(Timer::REAL);
        solution["solution"]["totalBest"] = bestSolution->scv;
        solution["solution"]["feasible"] = true;

        solution["solution"]["timeslots"] = arrayValue;
        for(int i = 0; i < (*bestSolution).data->n_of_events; i++) //{
            solution["solution"]["timeslots"].append(bestSolution->sln[i].first);

        solution["solution"]["rooms"] = arrayValue;
        for(int i = 0; i < (*bestSolution).data->n_of_events; i++) //{
            solution["solution"]["rooms"].append(bestSolution->sln[i].second);

        string solutionStr = writeString(swb, solution);
        (*os) << solutionStr << endl;
    }
    else {
        solution["solution"]["totalTime"] = timer.elapsedTime(Timer::REAL);
        solution["solution"]["totalBest"] = (bestSolution->computeHcv() * 1000000) + bestSolution->computeScv();
        solution["solution"]["feasible"] = false;

        string solutionStr = writeString(swb, solution);
        (*os) << solutionStr << endl;
    }
}

/*
 * Sets the bestScv and bestEvaluation
 * and prints the process's bestScv (if feasible) or bestEvaluation (if not feasible)
 */
void setCurrentCost(Solution *currentSolution, int tid, StreamWriterBuilder swb, Value logEntry) {//
    logEntry["logEntry"]["threadID"] = tid;
    int currentScv = currentSolution->scv;

    if(currentSolution->feasible) {
        if (currentScv != bestScv) {
            bestScv = currentScv;
            bestEvaluation = currentScv;
            logEntry["logEntry"]["best"] = bestScv;
            double time = getTime();
            logEntry["logEntry"]["time"] = time < 0 ? 0.0 : time;
            string logEntryStr = writeString(swb, logEntry);
            (*os) << logEntryStr << endl;
        }
    } else if(!currentSolution->feasible) {
        int currentEvaluation = (currentSolution->computeHcv() * 1000000 + currentSolution->computeScv()) ;
        if(currentEvaluation < bestEvaluation) {
            bestEvaluation = currentEvaluation;
            logEntry["logEntry"]["best"] = bestEvaluation;
            double time = getTime();
            logEntry["logEntry"]["time"] = time < 0 ? 0.0 : time;
            string logEntryStr = writeString(swb, logEntry);
            (*os) << logEntryStr << endl;
        }
    }
}

/*
 * Sets the global bestScv and bestEvaluation
 * and prints the global bestScv (if feasibel) or bestEvaluation (if not feasible)
 */
void setGlobalCost(Solution *currentSolution, StreamWriterBuilder swb, Value runEntry) {
    int currentScv = currentSolution->scv;
    if(currentSolution->feasible) {
        MPI_Allreduce(&currentScv, &bestScv, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        runEntry["runEntry"]["totalBest"] = bestScv;
        runEntry["runEntry"]["feasible"] = true;
        string runEntryStr = writeString(swb, runEntry);

        if (id == 0) {
            (*os) << runEntryStr << endl;
        }
    }
    else if(!currentSolution->feasible) {
        int currentEvaluation = (currentSolution->computeHcv() * 1000000 + currentSolution->computeScv());
        MPI_Allreduce(&currentEvaluation, &bestEvaluation, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        runEntry["runEntry"]["totalBest"] = bestEvaluation;
        runEntry["runEntry"]["feasible"] = false;
        string runEntryStr = writeString(swb, runEntry);

        if (id == 0) {
            (*os) << runEntryStr << endl;
        }
    }
}

/*
 * Serialize the Problem instance so it can be send as contiguous MPI message
 *
 * @param  outBuffer The out buffer that will hold the serialized problem
 */
void serializeProblem(void *outBuffer) {
    int sizeOf_studentNumbern = problem->n_of_events;
    int sizeOf_roomSize = problem->n_of_rooms;
    int sizeOf_student_events = problem->n_of_students*problem->n_of_events;
    int sizeOf_eventCorrelations = problem->n_of_events*problem->n_of_events;
    int sizeOf_room_features = problem->n_of_rooms*problem->n_of_features;
    int sizeOf_event_features = problem->n_of_events*problem->n_of_features;
    int sizeOf_possibleRooms = problem->n_of_events*problem->n_of_rooms;

    int position = 0;
    MPI_Pack(&problem->n_of_events, 1, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
    MPI_Pack(&problem->n_of_rooms, 1, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
    MPI_Pack(&problem->n_of_features, 1, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
    MPI_Pack(&problem->n_of_students, 1, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
    MPI_Pack(&problem->studentNumber.front(), sizeOf_studentNumbern, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
    MPI_Pack(&problem->roomSize.front(), sizeOf_roomSize, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
    MPI_Pack(*problem->student_events, sizeOf_student_events, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
    MPI_Pack(*problem->eventCorrelations, sizeOf_eventCorrelations, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
    MPI_Pack(*problem->room_features, sizeOf_room_features, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
    MPI_Pack(*problem->event_features, sizeOf_event_features, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
    MPI_Pack(*problem->possibleRooms, sizeOf_possibleRooms, MPI_INT, outBuffer, sizeOfProblem, &position, MPI_COMM_WORLD);
}

/*
 * Deserialize the Problem instance from a buffer
 *
 * @param  outBuffer The buffer that holds the data for the problem to be deserialized
 */
void deserializeProblem(void *bufferProblem) {
    int position = 0;
    MPI_Unpack(bufferProblem, sizeOfProblem, &position, &n_of_events, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(bufferProblem, sizeOfProblem, &position, &n_of_rooms, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(bufferProblem, sizeOfProblem, &position, &n_of_features, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(bufferProblem, sizeOfProblem, &position, &n_of_students, 1, MPI_INT, MPI_COMM_WORLD);

    problem = new Problem(n_of_events, n_of_rooms, n_of_features, n_of_students);

    MPI_Unpack(bufferProblem, sizeOfProblem, &position, &problem->studentNumber.front(), problem->n_of_events, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(bufferProblem, sizeOfProblem, &position, &problem->roomSize.front(), problem->n_of_rooms, MPI_INT, MPI_COMM_WORLD);

    MPI_Unpack(bufferProblem, sizeOfProblem, &position, *problem->student_events, problem->n_of_students*problem->n_of_events, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(bufferProblem, sizeOfProblem, &position, *problem->eventCorrelations, problem->n_of_events*problem->n_of_events, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(bufferProblem, sizeOfProblem, &position, *problem->room_features, problem->n_of_rooms*problem->n_of_features, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(bufferProblem, sizeOfProblem, &position, *problem->event_features, problem->n_of_events*problem->n_of_features, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(bufferProblem, sizeOfProblem, &position, *problem->possibleRooms, problem->n_of_events*problem->n_of_rooms, MPI_INT, MPI_COMM_WORLD);
}

/*
 * Serialize soluntions from the pop array so they can be send as contiguous MPI message
 *
 * @param  offset An offset that determines from where to begin serializing in the pop array
 * @param  numSolutions Number of sollutions to be serialized from the pop array
 * @param  outBuffer The out buffer that will hold the serialized sollutions
 */
void serializeSolutions(int offset, int numSolutions, void *outBuffer) {
    int position = 0;
    for(int i = offset; i < offset + numSolutions; i++) {
        int *slnBuffer = (int *)malloc(2 * pop[i]->data->n_of_events * sizeof(int));
        int k = -1;
        for (int j = 0; j < pop[i]->data->n_of_events; j++) {
            k = k + 1;
            slnBuffer[k] = pop[i]->sln[j].first;
            k = k + 1;
            slnBuffer[k] = pop[i]->sln[j].second;
        }
        MPI_Pack(slnBuffer, 2 * pop[i]->data->n_of_events, MPI_INT, outBuffer, numSolutions*sizeOfSolution, &position, MPI_COMM_WORLD);
        MPI_Pack(&pop[i]->feasible, 1, MPI_C_BOOL, outBuffer, numSolutions*sizeOfSolution, &position, MPI_COMM_WORLD);
        MPI_Pack(&pop[i]->scv, 1, MPI_INT, outBuffer, numSolutions*sizeOfSolution, &position, MPI_COMM_WORLD);
        MPI_Pack(&pop[i]->hcv, 1, MPI_INT, outBuffer, numSolutions*sizeOfSolution, &position, MPI_COMM_WORLD);
        MPI_Pack(&pop[i]->penalty, 1, MPI_INT, outBuffer, numSolutions*sizeOfSolution, &position, MPI_COMM_WORLD);
    }
}

/*
 * Deserialize soluntions from a buffer so they can be added to the pop array
 *
 * @param  offset An offset that determines from where to begin deserializing to the pop array
 * @param  numSolutions Number of sollutions to be deserialized
 * @param  sollutionsBuffer The buffer from where to deserialize the solutions
 */
void deserializeSolution(int offset, int numSolutions, void *sollutionsBuffer) {
    int position = 0;
    for(int i = popSize - numSolutions - offset; i < popSize - offset; i++) {
        pop[i] = new Solution(problem, rnd);
        int *slnBuffer = (int *)malloc(2 * pop[i]->data->n_of_events * sizeof(int));
        MPI_Unpack(sollutionsBuffer, sizeOfSolution, &position, slnBuffer, 2 * pop[i]->data->n_of_events, MPI_INT, MPI_COMM_WORLD);
        pair<int,int> initPair;
        int k = -1;
        for (int j = 0; j < pop[i]->data->n_of_events; j++) {
            k = k + 1;
            initPair.first = slnBuffer[k];
            k = k + 1;
            initPair.second = slnBuffer[k];
            pop[i]->sln[j] = initPair;
        }
        MPI_Unpack(sollutionsBuffer, sizeOfSolution, &position, &pop[i]->feasible, 1, MPI_C_BOOL, MPI_COMM_WORLD);
        MPI_Unpack(sollutionsBuffer, sizeOfSolution, &position, &pop[i]->scv, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(sollutionsBuffer, sizeOfSolution, &position, &pop[i]->hcv, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(sollutionsBuffer, sizeOfSolution, &position, &pop[i]->penalty, 1, MPI_INT, MPI_COMM_WORLD);
        for (int j = 0; j < pop[i]->data->n_of_events; j++ ){
            int t = pop[i]->sln[j].first;
            pop[i]->timeslot_events[t].push_back(j);
        }
    }
}

int main( int argc, char** argv) {

    // Initialize MPI execution
    rc = MPI_Init(&argc, &argv);
    if (rc != 0) {
        (*os) << "ERROR MPI" << endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (id == 0) {
        // The master (server) process

        Control control(argc, argv);
        problemType = control.getProblemType();

        if (problemType == 1){
            maxSteps = 200;
        }
        else if (problemType == 2) {
            maxSteps = 1000;
        }
        else{
            maxSteps = 2000;
        }

        problem = new Problem(control.getInputStream());

        seed = (unsigned) control.getSeed();
        rnd = new Random(seed);

        t = control.getThreadsNum();
        MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // broadcast the maxSteos for the local search, determined from the problem type
        MPI_Bcast(&maxSteps,  1, MPI_INT, 0, MPI_COMM_WORLD);

        for(int i = 1; i < p; i++) {
            // calculate every process's seed
            int proc_seed = abs(control.getSeed() + i*(control.getSeed()/10));
            // the master process send every process its seed
            MPI_Send(&proc_seed, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        sizeOfProblem = getProblemSizeinB(problem);
        bufferProblem = malloc((size_t)sizeOfProblem);

        // broadcast the problem size - needed for the migration
        MPI_Bcast(&sizeOfProblem,  1, MPI_INT, 0, MPI_COMM_WORLD);

        serializeProblem(bufferProblem);

        // broadcast the problem instance
        MPI_Bcast(bufferProblem, sizeOfProblem, MPI_PACKED, 0, MPI_COMM_WORLD);

        // generate an initial population
        for(int i=0; i < popSize; i++){
            pop[i] = new Solution(problem, rnd);
            pop[i]->RandomInitialSolution();
            pop[i]->localSearch(maxSteps);
            pop[i]->computePenalty();
        }

        sizeOfSolution = getSolutionSizeInB(problem);
        bufferSolution = malloc((size_t)(popSize*sizeOfSolution));

        // broadcast the solution size - needed fro the migration
        MPI_Bcast(&sizeOfSolution,  1, MPI_INT, 0, MPI_COMM_WORLD);

        serializeSolutions(0, popSize, bufferSolution);
        // broadcast the initial population
        MPI_Bcast(bufferSolution, (size_t)popSize*sizeOfSolution, MPI_PACKED, 0, MPI_COMM_WORLD);
    } else {
        // The other processes in the communicator (except the master process)

        MPI_Bcast(&t,  1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Bcast(&maxSteps,  1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Status status;
        MPI_Recv(&seed, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        rnd = new Random(seed);

        MPI_Bcast(&sizeOfProblem,  1, MPI_INT, 0, MPI_COMM_WORLD);
        bufferProblem = malloc((size_t)sizeOfProblem);
        MPI_Bcast(bufferProblem, sizeOfProblem, MPI_PACKED, 0, MPI_COMM_WORLD);
        deserializeProblem(bufferProblem);

        MPI_Bcast(&sizeOfSolution,  1, MPI_INT, 0, MPI_COMM_WORLD);
        bufferSolution = malloc((size_t)popSize*sizeOfSolution);
        MPI_Bcast(bufferSolution, (size_t)popSize*sizeOfSolution, MPI_PACKED, 0, MPI_COMM_WORLD);
        deserializeSolution(0, popSize, bufferSolution);
    }

    // Continue execution from all processes (the master and the other processes play a client role)

    StreamWriterBuilder swb;
    swb.settings_["indentation"] = "";
    Value runEntry;
    Value logEntry;
    Value solution;
    solution["solution"]["procID"] = id;

    beginTry();

    int generation;
    int snd = (id + 1) % p; // determine the id of the process to send solutions to
    int rcv = (p + id - 1) % p; // determine the id of the process from which to receive solutions
    int numMigrants = 1; // number of migrants for each of the two direction migrations

    void *sendBuffer;
    void *recvBuffer;

    int numberMigrationPeriods = 0;

    omp_set_num_threads(t);

#pragma omp parallel private(generation, numberMigrationPeriods, tid, logEntry)//t, child, copyParent1, copyParent2, swb
    {
        Solution *child;
        Solution *copyParent1;
        Solution *copyParent2;

        numberMigrationPeriods = 0;

        tid = omp_get_thread_num();

        StreamWriterBuilder swb;
        swb.settings_["indentation"] = "";
        logEntry["logEntry"]["procID"] = id;
        setCurrentCost(pop[0], tid, swb, logEntry);

        if (tid == 0) {
            sendBuffer = malloc((size_t)(sizeOfSolution * numMigrants));
            recvBuffer = malloc((size_t)(sizeOfSolution * numMigrants));
        }

        for (generation = tid; generation <= 2000; generation += t) {
            numberMigrationPeriods++;

            // migration on every 10-th part of every thread's generations (with an offset of 50 iterations from the beginning)
            if (numberMigrationPeriods > 0 && numberMigrationPeriods % 100 == 50) {
                // the migration is performed only from the master thread
                if (tid == 0) {
#pragma omp critical
                    {
                        // MPI_Barrier for beginning of synchronization of the master threads of all MPI processes
                        MPI_Barrier(MPI_COMM_WORLD);

                        serializeSolutions(0, numMigrants, sendBuffer);

                        MPI_Status status;
                        MPI_Sendrecv(sendBuffer, (sizeOfSolution * numMigrants), MPI_PACKED, snd, SEND_RECV_TAG,
                                     recvBuffer, (sizeOfSolution * numMigrants), MPI_PACKED, rcv, SEND_RECV_TAG, MPI_COMM_WORLD, &status);

                        deserializeSolution(0, numMigrants, recvBuffer);

                        serializeSolutions(numMigrants, numMigrants, sendBuffer);

                        MPI_Sendrecv(sendBuffer, (sizeOfSolution * numMigrants), MPI_PACKED, rcv, SEND_RECV_TAG_2,
                                     recvBuffer, (sizeOfSolution * numMigrants), MPI_PACKED, snd, SEND_RECV_TAG_2, MPI_COMM_WORLD, &status);

                        deserializeSolution(numMigrants, numMigrants, recvBuffer);

                        // MPI_Barrier for end of synchronization of the master threads of all MPI processes
                        MPI_Barrier(MPI_COMM_WORLD);
                    }
                }
            }

            child = new Solution(problem,rnd);
            child->RandomInitialSolution();
            copyParent1 = new Solution(problem,rnd);
            copyParent1->RandomInitialSolution();
            copyParent2 = new Solution(problem,rnd);
            copyParent2->RandomInitialSolution();

            // select parents
            Solution *parent1 = selection5(pop, popSize);
            Solution *parent2 = selection5(pop, popSize);

            // copy the parents to local thread vars so that the omp threads do not interfere with one another
#pragma omp critical
            {
                copyParent1->copy(parent1);
                copyParent2->copy(parent2);
            }

            // generate child
            if(rnd->next() < 0.8)
                child->crossover(copyParent1,copyParent2);
            else {
                child = copyParent1;
            }

            // do some mutation
            if(rnd->next() < 0.5) {
                child->mutation();
            }

            //apply local search to offspring
            child->localSearch(maxSteps);

            //evaluate the offspring
            child->computePenalty();

            // critical section on the shared popilation array - replace worst member of the population with offspring
#pragma omp critical
            {
                pop[popSize - 1]->copy(child);
                sort(pop, pop + popSize, compareSolution);
                setCurrentCost(pop[0], tid, swb, logEntry);
            }
            delete child;
        }
    }

    setGlobalCost(pop[0], swb, runEntry);

    endTry(pop[0], solution);

    // remember to delete the population
    for(int i=0; i < popSize; i++){
        delete pop[i];
    }

    delete problem;
    delete rnd;

    elapsed_time += MPI_Wtime();
    if (!id) {
        runEntry["runEntry"]["procsNum"] = p;
        runEntry["runEntry"]["threadsNum"] = t;
        runEntry["runEntry"]["totalTime"] = elapsed_time;
        string runEntryStr = writeString(swb, runEntry);
        (*os) << runEntryStr << endl;
    }
    //Finalize MPI execution
    MPI_Finalize();
    return 0;
}
