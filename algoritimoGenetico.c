// file: geneticAlgorithm1.c
// by: John LeFlohic
// use: simulates a species' convergence to a model
// made: February 24, 1999

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include "mpi.h"

#define NUMBER_ORGANISMS 10000
#define NUMBER_GENES     50
#define ALLELES          4
#define MUTATION_RATE    0.001

#define MAXIMUM_FITNESS NUMBER_GENES
#define FALSE 0
#define TRUE  1


// global variables
char **currentGeneration, **nextGeneration;
char modelOrganism[NUMBER_GENES];
int *organismsFitnesses;
int totalOfFitnesses;


// function declarations
void    AllocateMemory(void);
int     DoOneRun(void);
void    InitializeOrganisms(void);
int     EvaluateOrganisms(void);
void    ProduceNextGeneration(void);
int     SelectOneOrganism(void);
void    CreateModel(void);

// functions
int main()
{

    int ret, rank, size, i, tag, finalGeneration;
    double tempo;
    /*
    MPI_Status status;


    ret = MPI_Init(&argc, &argv);
    ret = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ret = MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (ret < 0)
    {
        printf ("Nao foi possivel inicializar o processo MPI!\n");
        exit();
    }
    else
    {
        tempo=MPI_Wtime();

    */
        //if(rank == 0){
            CreateModel();
        //}

        //for(i = 1; i < size; i++){
          //  if(rank == i){
                AllocateMemory();
                finalGeneration = DoOneRun();
                printf("\nThe final generation was: %d\n", finalGeneration);
            //}
        //}

        /*
        tempo=MPI_Wtime()-tempo;
        printf("Tempo decorrido no processo %d = %f\n", rank, tempo);
        MPI_Finalize();
    }
    */

}

//cria o organismo modelo
void CreateModel(void){
    int gene;
    for(gene=0; gene<NUMBER_GENES; ++gene)
    {
        modelOrganism[gene] = rand()%ALLELES;
    }
}

void AllocateMemory(void)
{
    printf("\nalocando memoria");
    int organism;

    currentGeneration =  (char**)malloc(sizeof(char*) * NUMBER_ORGANISMS);
    nextGeneration =     (char**)malloc(sizeof(char*) * NUMBER_ORGANISMS);
    //modelOrganism =      (char*)malloc(sizeof(char) * NUMBER_GENES);
    organismsFitnesses = (int*)malloc(sizeof(int) * NUMBER_ORGANISMS);

    for(organism=0; organism<NUMBER_ORGANISMS; ++organism)
    {
        currentGeneration[organism] =   (char*)malloc(sizeof(char) * NUMBER_GENES);
        nextGeneration[organism] =      (char*)malloc(sizeof(char) * NUMBER_GENES);
    }
}


int DoOneRun(void)
{
    int generations = 1;
    int perfectGeneration = FALSE;

    InitializeOrganisms();

    while(TRUE)
    {
        perfectGeneration = EvaluateOrganisms();
        if( perfectGeneration==TRUE )
            return generations;
        ProduceNextGeneration();
        ++generations;
    }
}


void InitializeOrganisms(void)
{
    printf("\niniciando orgranismos");
    int organism;
    int gene;

    // initialize the normal organisms
    for(organism=0; organism<NUMBER_ORGANISMS; ++organism)
    {
        for(gene=0; gene<NUMBER_GENES; ++gene)
        {
            currentGeneration[organism][gene] = rand()%ALLELES;
        }
    }
}


int EvaluateOrganisms(void)
{
    //printf("\navaliando organismos");
    int organism;
    int gene;
    int currentOrganismsFitnessTally;

    totalOfFitnesses = 0;

    for(organism=0; organism<NUMBER_ORGANISMS; ++organism)
    {
        currentOrganismsFitnessTally = 0;

        // tally up the current organism's fitness
        for(gene=0; gene<NUMBER_GENES; ++gene)
        {
            if( currentGeneration[organism][gene]
                    == modelOrganism[gene] )
            {
                ++currentOrganismsFitnessTally;
            }
        }

        // save the tally in the fitnesses data structure
        // and add its fitness to the generation's total
        organismsFitnesses[organism] =
            currentOrganismsFitnessTally;
        totalOfFitnesses += currentOrganismsFitnessTally;

        // check if we have a perfect generation
        if( currentOrganismsFitnessTally == MAXIMUM_FITNESS )
        {
            return TRUE;
        }
    }
    return FALSE;

}


void ProduceNextGeneration(void)
{
    //printf("\nproduzindo nova geracao");
    int organism;
    int gene;
    int parentOne;
    int parentTwo;
    int crossoverPoint;
    int mutateThisGene;

    // fill the nextGeneration data structure with the
    // children
    for(organism=0; organism<NUMBER_ORGANISMS; ++organism)
    {
        parentOne = SelectOneOrganism();
        parentTwo = SelectOneOrganism();
        crossoverPoint  =  rand() % NUMBER_GENES;

        for(gene=0; gene<NUMBER_GENES; ++gene)
        {

            // copy over a single gene
            mutateThisGene = rand() % (int)(1.0 / MUTATION_RATE);
            if(mutateThisGene == 0)
            {

                // we decided to make this gene a mutation
                nextGeneration[organism][gene] = rand() % ALLELES;
            }
            else
            {
                // we decided to copy this gene from a parent
                if (gene < crossoverPoint)
                {
                    nextGeneration[organism][gene] =
                        currentGeneration[parentOne][gene];
                }
                else
                {
                    nextGeneration[organism][gene] =
                        currentGeneration[parentTwo][gene];
                }
            }
        }
    }

    // copy the children in nextGeneration into
    // currentGeneration
    for(organism=0; organism<NUMBER_ORGANISMS; ++organism)
    {
        for(gene=0; gene<NUMBER_GENES; ++gene)
        {
            currentGeneration[organism][gene] =
                nextGeneration[organism][gene];
        }
    }
}


int SelectOneOrganism(void)
{
    //printf("\nseleciona um orgranismo");
    int organism;
    int runningTotal;
    int randomSelectPoint;

    runningTotal = 0;
    randomSelectPoint = rand() % (totalOfFitnesses + 1);

    for(organism=0; organism<NUMBER_ORGANISMS; ++organism)
    {
        runningTotal += organismsFitnesses[organism];
        if(runningTotal >= randomSelectPoint)
            return organism;
    }
}
