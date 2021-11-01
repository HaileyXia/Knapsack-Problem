//
//  main.c
//  mknapsack
//
//  Created by Bai on 14/03/2020.
//  Copyright © 2019 UNNC. All rights reserved.
//  PSO algorithm for MKP 
// 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

/* global parameters */
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
int NUM_OF_RUNS = 1; //number of runs for a single problem
int MAX_TIME = 300;  //max amount of time permited (in sec)
int num_of_problems; //number of problems

/* parameters for PSO*/   
clock_t time_start, time_finish;
static int SWARM_SIZE = 1000; 
float cognitive_learning = 7; 
float social_learning = 0.9;  
float inertia = 1; 
int MAX_NUM_OF_ITER = 10000; //max number of generations
int fitnessP = 5000; //parameter for the fitness func
float v_max = 4;     //max velocity

struct move_struct{
    int item1;
    int item2;
};

struct solution_struct best_sln;  //global best solution

//return a random number between 0 and 1 
float rand_01()
{
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    //printf("rand01=%f\n", number);
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max)
{
    int div = max-min+1;
    int val =rand() % div + min;
    //printf("rand_range= %d \n", val);
    return val;
}

struct item_struct{
    int dim; //no. of dimensions
    int* size; //volume of item in all dimensions
    int p; //profit
    double ratio; //profit/volume
    int indx;  //index of the item
};

struct problem_struct{
    int n; //number of items
    int dim; //number of dimensions
    struct item_struct* items;
    int* capacities;  //knapsack capacities
    int optimal;       //optimal solution
};

//free memory
void free_problem(struct problem_struct* prob)
{
    if(prob!=NULL)
    {
        if(prob->capacities !=NULL) free(prob->capacities);
        if(prob->items!=NULL)
        {
            for(int j=0; j<prob->n; j++)
            {
                if(prob->items[j].size != NULL)
                    free(prob->items[j].size);
            }
            free(prob->items);
        }
        free(prob);
    }
}

//initialize problem, allocate memory
void init_problem(int n, int dim, struct problem_struct** my_prob)
{
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    new_prob->n=n; new_prob->dim=dim;
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    for(int j=0; j<n; j++)
        new_prob->items[j].size= malloc(sizeof(int)*dim);
    new_prob->capacities = malloc(sizeof(int)*dim);
    *my_prob = new_prob;
}


//example to create problem instances, actual date should come from file
struct problem_struct** load_problems(char* data_file)
{
    int i,j,k;
    //int num_of_probs;
    FILE* pfile = fopen(data_file, "r");
    if(pfile==NULL)
        {printf("Data file %s does not exist. Please check!\n", data_file); exit(2); }
    fscanf (pfile, "%d", &num_of_problems);
 
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct*)*num_of_problems);
    for(k=0; k<num_of_problems; k++)
    {
        int n, dim, obj_opt;
        fscanf (pfile, "%d", &n); //item num
        fscanf (pfile, "%d", &dim);  //dim num
        fscanf (pfile, "%d", &obj_opt); 

        init_problem(n, dim, &my_problems[k]);  //allocate data memory

        my_problems[k]->optimal = obj_opt;   //

        for(j=0; j<n; j++)
        {
            my_problems[k]->items[j].dim=dim;
            my_problems[k]->items[j].indx=j;
            fscanf(pfile, "%d", &(my_problems[k]->items[j].p)); //read profit data
            //printf("item[j].p=%d\n",my_problems[k]->items[j].p);
        }
        for(i=0; i<dim; i++)
        {
            for(j=0; j<n; j++)
            {
                fscanf(pfile, "%d", &(my_problems[k]->items[j].size[i])); //read size data
                //printf("my_problems[%i]->items[%i].size[%i]=%d\n",k,j,i,my_problems[k]->items[j].size[i]);
            }
        }
        for(i=0; i<dim; i++){
            fscanf(pfile, "%d", &(my_problems[k]->capacities[i]));
            //printf("capacities[i]=%d\n",my_problems[k]->capacities[i] );
        }
    }
    fclose(pfile); //close file
    return my_problems;
}

struct solution_struct{
    struct problem_struct* prob; //maintain a shallow copy of the problem data
    float objective;   //the price of the current solution
    int feasibility; //indicate the feasibility of the solution
    int* x; // solution encoding vector
    int* cap_left; //capacity left in all dimensions
    float* velocity; //velocity vector
    float fitness; //fitness value for each particle(solution)
};

//free solution memory
void free_solution(struct solution_struct* sln)
{
    if(sln!=NULL)
    {
        free(sln->x);
        free(sln->cap_left);
        free(sln->velocity);
        sln->fitness=0;
        sln->objective=0;
        sln->prob=NULL;
        sln->feasibility=false;
    }
}

//copy a solution from another solution
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln)
{
    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {
        dest_sln = malloc(sizeof(struct solution_struct));
    }
    else{
        free(dest_sln->cap_left);
        free(dest_sln->x);
    }
    int n = source_sln->prob->n;
    int m =source_sln->prob->dim;
    dest_sln->x = malloc(sizeof(int)*n); //allocate memory to target solution
    dest_sln->cap_left=malloc(sizeof(int)*m);
    for(int i=0; i<m; i++)
        dest_sln->cap_left[i]= source_sln->cap_left[i];
    for(int j=0; j<n; j++)
        dest_sln->x[j] = source_sln->x[j];
    dest_sln->prob= source_sln->prob;
    dest_sln->feasibility=source_sln->feasibility;
    dest_sln->objective=source_sln->objective;
    return true;
}

//evaluate the solution, check valid or not
void evaluate_solution(struct solution_struct* sln) 
{
    //evaluate the feasibility and objective of the solution
    sln->objective =0; sln->feasibility = 1;
    struct item_struct* items_p = sln->prob->items;
    //printf("a\n");
    for(int i=0; i< items_p->dim; i++) //check in all dimension
    {
        sln->cap_left[i]=sln->prob->capacities[i];  
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->cap_left[i] -= items_p[j].size[i]*sln->x[j];
            if(sln->cap_left[i]<0) {
                sln->feasibility = -1*i; //if exceeding capacity, obj=0
                return;
            }
        }
    }
    if(sln->feasibility>0) //if solution is valid, calculate the objective value
    {
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->objective += sln->x[j] * items_p[j].p;
        }
    }
}

//output a given solution to a file
void output_solution(struct solution_struct* sln, char* out_file)
{
    if(out_file !=NULL){
        FILE* pfile = fopen(out_file, "a"); //append solution data
        fprintf(pfile, "%i\n", (int)sln->objective);
        for(int i=0; i<sln->prob->n; i++)
        {
            //printf("%d\n",sln->x[i]);
            fprintf(pfile, "%i ", sln->x[i]);
        }
        fprintf(pfile, "\n");
        /*for(int j=0; j<sln->prob->n; j++)
            fprintf(pfile, "%i ", sln->prob->items[j].p);
        fprintf(pfile, "\n");*/
        fclose(pfile);
    }
    else
        printf("sln.feas=%d, sln.obj=%f\n", sln->feasibility, sln->objective);
}


//check the  feasiblity and obj values of solutons from solution_file.
//return 0 is all correct or the index of the first infeasible problem [1, num_of_problems].
int check_solutions(struct problem_struct** my_problems, char* solution_file)
{
    FILE * pfile= fopen(solution_file, "r");
    if(pfile==NULL)
    {
        printf("Solution file %s does not exist. Please check!\n", solution_file);
        exit(2);
    }
    float val_obj;
    int val;
    fscanf (pfile, "%i", &val);
    if(val != num_of_problems)
    {
        printf("The stated number of solutions does not match the number of problems.\n");
        exit(3);
    }
    struct solution_struct temp_sln;
    int count=0, k=0;
    int n, dim;
    while(fscanf (pfile, "%f", &val_obj)!=EOF && k<num_of_problems)
    {
        //val_obj = val;
        n= my_problems[k]->n;  dim= my_problems[k]->dim;
        temp_sln.x = malloc(sizeof(int)*n);
        temp_sln.cap_left=malloc(sizeof(int)*dim);
        temp_sln.prob = my_problems[k];
        while(fscanf (pfile, "%i", &val)!=EOF)
        {
            if(val<0 || val>1) {fclose(pfile);  return k+1;} //illeagal values
            temp_sln.x[count] = val;
            count++;
            if(count==n)
            {
                evaluate_solution(&temp_sln);
                if(!temp_sln.feasibility || fabs(temp_sln.objective - val_obj)>0.01)
                {
                    fclose(pfile);
                    //printf("feasb=%i, obj= %f, val=%i\n",temp_sln.feasibility, temp_sln.objective, val_obj);
                    //output_solution(&temp_sln, "my_debug.txt");
                    return k+1;  //infeasible soltuion or wrong obj
                }
                else{
                    break;
                }
            }
        }
        count=0; k++;
        
        free(temp_sln.x); free(temp_sln.cap_left);
    }
    fclose(pfile);
    return 0;
}


// add item back to package
void add_item(struct solution_struct* pop, struct item_struct* sort_item){
    int item_num = pop->prob->n;
    int max_index;
    while (pop->feasibility == 1){ //continue adding until package is full
        for (int i = 0; i < item_num; i++){ //initialize adding item
            if (pop->x[i] == 0){
                max_index = i;
                break;
            }
        }
        
        for(int j = 0; j < item_num; j++){ //get available item with largest ratio
            if(pop->x[sort_item[j].indx] == 0){
                max_index = sort_item[j].indx;
                break;
            }
        }
        pop->x[max_index] = 1; //add item with largest ratio
        evaluate_solution(pop);   
    }
    pop->x[max_index] = 0;  //last one is over capacity, remove it
    evaluate_solution(pop);  
}


//modify the solutions that violate the capacity constraints, based on ratio
void feasibility_repair(struct solution_struct* pop, struct item_struct* sort_item)
{
    int item_num = pop->prob->n;
    while(pop->feasibility != 1){
        int min_index;
        for(int j = 0; j < item_num; j++){
            if(pop->x[j] == 1){
                min_index = j;
                break;
            }
        }
        //initialize the min_index
        for(int j = item_num-1; j >= 0; j--){
            if(pop->x[sort_item[j].indx] == 1){
                min_index = sort_item[j].indx;
                break;
            }
        }
        pop->x[min_index] = 0; //remove the item with minimal ratio
        evaluate_solution(pop); 
    }
}

//update global best solution from sln
void update_best_solution(struct solution_struct* sln)
{
    if(best_sln.objective < sln->objective)
    copy_solution(&best_sln, sln);
}

// calculate the fitness value of each particle(solution)
int fitnessFunction(struct solution_struct* sln){ 
    int total_poslin = 0;
    int profit = 0;
    sln->fitness = 0;
    int poslin = 0;
    struct item_struct* item_i = sln->prob->items;
    for (int i = 0; i < item_i->dim; i++) //cap_left in each dimension
    {
        if (sln->cap_left[i]>=0)  //valid solution
        {
            poslin =0;
        }
        else{   //invalid solution
            poslin = -(sln->cap_left[i]);
        }
        total_poslin = total_poslin + poslin;
    }
    for (int i = 0; i < sln->prob->n; i++) //calculate the total profit
    {
        profit = profit + sln->x[i] * item_i[i].p;
    }
    sln->fitness = profit - fitnessP * total_poslin; //get the fitness value
    return sln->fitness;
}

//intialise the swarm with random solutions 
struct solution_struct* init_swarm(struct problem_struct* prob, struct solution_struct* pop){   
    for(int i = 0; i < SWARM_SIZE; i++){
        pop[i].prob = prob;
        pop[i].x = malloc(sizeof(int)*prob->n);
        pop[i].velocity = malloc(sizeof(int)*prob->n);
        pop[i].cap_left = malloc(sizeof(int)*prob->dim);    
        //initialize all items as unpacked and initialize velocity
        for(int j = 0; j < prob->n; j++){
            pop[i].x[j] = 0;
            pop[i].velocity[j] = rand_int(-1*v_max,v_max);
        }
        //initialize capacities in all dimensions
        for(int k = 0; k < prob->dim; k++){
            pop[i].cap_left[k]=prob->capacities[k];
        }

        //randomly initialize x      
        while(1){
            int index = rand_int(0, prob->n-1);
            pop[i].x[index] = 1;
            bool space_left = true;
            for(int k = 0; k < prob->dim; k++){
                pop[i].cap_left[k] -= prob->items[index].size[k]; //check in each dimension
                if(pop[i].cap_left[k] < 0){
                    space_left = false;
                }
            }
            //unpack the last item
            if(!space_left){   
                pop[i].x[index] = 0;
                for(int k = 0; k < prob->dim; k++){
                    pop[i].cap_left[k] += prob->items[index].size[k];
                    //printf("cap_left: %d\n",pop[i].cap_left[k]);
                }
                break;
            }
        }
        //printf("%f\n",pop[i].objective);
        evaluate_solution(&pop[i]);  //check the solution is valid, get objective value and feasibility value
        //printf("randomInit_sln obj=%f\tfeasiblity = %d, \n", pop[i].objective, pop[i].feasibility);
    }
    return pop;
}

//compare function
int cmpfunc1(const void* a, const void* b){
    const struct item_struct* item1 = a;
    const struct item_struct* item2 = b;
    if(item1->ratio>item2->ratio) return -1;
    if(item1->ratio<item2->ratio) return 1;
    return 0;
    }
int cmpfunc2 (const void * a, const void * b) {
        const struct item_struct* item1 = a;
        const struct item_struct* item2 = b;
        if(item1->indx>item2->indx) return 1;
        if(item1->indx<item2->indx) return -1;
        return 0;
    }

int cmpfunc_sln (const void * a, const void * b) {
    const struct solution_struct* sln1 = a;
    const struct solution_struct* sln2 = b;
    if(sln1->objective > sln2 ->objective) return -1;
    if(sln1->objective < sln2 ->objective) return 1;
    return 0;
}

//a greedy heuristic solution based on profit/volume ratio
struct solution_struct* greedy_heuristic(struct problem_struct* prob)
{
    //sort the items based on ratio, desending order
    for(int i=0; i<prob->n;i++){
        double avg_size=0;
        struct item_struct* item_i = &prob->items[i];
        for(int d=0; d< prob->dim; d++){
            avg_size += (double)item_i->size[d]/prob->capacities[d];
        }
        item_i->ratio = item_i->p/avg_size;
    }
    qsort(prob->items, prob->n, sizeof(struct item_struct), cmpfunc1);
    //allocate memory
    struct solution_struct* init_sln = malloc(sizeof(struct solution_struct));
    init_sln->prob=prob;    init_sln->objective =0;
    init_sln->x = malloc(sizeof(int)*prob->n);
    init_sln->velocity = malloc(sizeof(int)*prob->n);
    init_sln->cap_left = malloc(sizeof(int)*prob->dim);
    int* cap = malloc(sizeof(int)*prob->dim);
    int i=0, d=0;
    for(d=0; d<prob->dim; d++) cap[d]=0; //aggregated volume
    for(i=0; i<prob->n; i++)
    {
        init_sln->velocity[i] = rand_int(-1*v_max,v_max);
        struct item_struct* item_i = &prob->items[i];
        //printf("item[%d].ratio = %.3f. profit %d\n",item_i->indx,prob->items[i].ratio, item_i->p);
        for(d=0; d<prob->dim; d++){
            //printf("%d\t",item_i->size[d]);
            if(cap[d] + item_i->size[d] > prob->capacities[d])
                break; //infeasible to pack this item, try next
        }
        if(d>=prob->dim){ //add item into package
            init_sln->x[item_i->indx] = 1;
            init_sln->objective += item_i->p;
            for(d=0; d<prob->dim; d++){
                cap[d] += item_i->size[d];
            }
            //printf("packing item %d\n", item_i->indx);
        }
        else{
            init_sln->x[item_i->indx] =0;
        }
    }
    for(d=0; d<prob->dim; d++){ //calculate cap left
        init_sln->cap_left[d] = prob->capacities[d]- cap[d];
    }
    free(cap);
    //restore item original order by sorting by index.
    qsort(prob->items, prob->n, sizeof(struct item_struct), cmpfunc2);
    
    evaluate_solution(init_sln); //check objective and feasibility of the solution
    //output_solution(init_sln, "greedy_sln.txt");
    //printf("Init_sln obj=\t%f\tfeasiblity = %d.\n", init_sln->objective, init_sln->feasibility);
    return init_sln;
}

// main PSO algorithm
int particleSwarm(struct problem_struct* prob){
    time_start = clock(); 
    double time_spent=0;
    int iter =0;        //number of interation
    struct solution_struct* curt_sln = greedy_heuristic(prob); // greedy initialization
    struct solution_struct swarm_pop[SWARM_SIZE];    
    struct solution_struct* pop = init_swarm(prob, swarm_pop); //random initialization
    struct solution_struct local_best[SWARM_SIZE];
    struct item_struct* items_particle = prob->items;
    struct item_struct sort_item[prob->n];
    
    // sort the item based on the profit/volume
    for(int i=0; i<prob->n;i++){
        double sort_size=0;
        struct item_struct* item_j = &prob->items[i];
        for(int d=0; d< prob->dim; d++){
            sort_size += (double)item_j->size[d]/prob->capacities[d];
        }
        item_j->ratio = item_j->p/sort_size;
    }
    qsort(prob->items, prob->n, sizeof(struct item_struct), cmpfunc1);

    for (int i = 0; i < prob->n; i++) 
    {
        sort_item[i] = prob->items[i]; //put sorted items into new list
    }
    //restore item original order by sorting by index.
    qsort(prob->items, prob->n, sizeof(struct item_struct), cmpfunc2);
    
    //allocate memory for the local best solution
    for(int i = 0; i < SWARM_SIZE; i++){
        local_best[i].x = malloc(sizeof(int)*pop->prob->n);
        local_best[i].cap_left = malloc(sizeof(int)*pop->prob->dim); 
    }

    copy_solution(&local_best[0],curt_sln);
    copy_solution(&pop[0],curt_sln);
    
    for (int i = 1; i < SWARM_SIZE; i++) //initialize the local optimal solutions
    {
        copy_solution(&local_best[i],&pop[i]);
    }
    
    best_sln.x = malloc(sizeof(int)*pop->prob->n); //initialize global optimal solution
    best_sln.cap_left = malloc(sizeof(int)*pop->prob->dim); 
    copy_solution(&best_sln,curt_sln);

    while(time_spent < MAX_TIME){   
        for (int i = 0; i < SWARM_SIZE; i++) 
        {
            if (pop[i].feasibility != 1)
            {
                feasibility_repair(&pop[i],sort_item); //repair all the invalid solutions
            }
            add_item(&pop[i],sort_item);  //try to add items into package if possible
        }
        
        for (int i = 0; i < SWARM_SIZE; i++) 
        {
            if (pop[i].objective>=local_best[i].objective) //update the local optimal solution
            {
                copy_solution(&local_best[i],&pop[i]); 
            }
            if (pop[i].objective>=best_sln.objective) //update the global optimal solution
            {
                copy_solution(&best_sln,&pop[i]); 
            }
            /*if (fitnessFunction(&pop[i])>=fitnessFunction(&local_best[i])) //call the fitness function
            {
                copy_solution(&local_best[i],&pop[i]);  
            }
            if (fitnessFunction(&pop[i])>=fitnessFunction(&best_sln))
            {
                copy_solution(&best_sln,&pop[i]);
            }*/

            for (int k = 0; k < prob->n; k++) //update the velocity and position
            {
                pop[i].velocity[k] = inertia*pop[i].velocity[k]+rand_01()*cognitive_learning*(local_best[i].x[k]-pop[i].x[k])+rand_01()*social_learning*(best_sln.x[k]-pop[i].x[k]); //这是在和上一个iteration里的比
                if (pop[i].velocity[k]>=v_max)
                {
                    pop[i].velocity[k] = v_max;
                }
                else if(pop[i].velocity[k]<=(-1*v_max)){
                    pop[i].velocity[k] = -1*v_max;
                }
                
                double sigma = rand_01();
                double solution_decode = 1.0/(1+exp(-pop[i].velocity[k]));
                
                if (solution_decode>=sigma)
                {
                    pop[i].x[k] = 1;
                }
                else{
                    pop[i].x[k] = 0;
                }
            }
            evaluate_solution(&pop[i]);  //evaluate the solution, get the objective value and feasibility 
        }
        iter++;
        time_finish=clock();
        time_spent = (double)(time_finish-time_start)/CLOCKS_PER_SEC;
    }
    update_best_solution(&best_sln); //update the best solution
    //printf("gap: %f  iter: %d  time: %f\n", (pop->prob->optimal-best_sln.objective)/pop->prob->optimal, iter, time_spent); 
    //printf("best: %f  iter: %d  time: %f\n", best_sln.objective, iter, time_spent);
    free_solution(curt_sln); //free memory
    free_solution(pop);
    return 0;
}

int main(int argc, const char * argv[]) {
    printf("Starting the run! \n");
    char data_file[50]={"somefile"}, out_file[50]={}, solution_file[50]={};  //max 50 problem instances per run
    if(argc<3)
    {
        printf("Insufficient arguments. Please use the following options:\n   -s data_file (compulsory)\n   -o out_file (default my_solutions.txt)\n   -c solution_file_to_check\n   -t max_time (in sec)\n");
        return 1;
    }
    else if(argc>9)
    {
        printf("Too many arguments.\n");
        return 2;
    }
    else
    {
        for(int i=1; i<argc; i=i+2)
        {
            if(strcmp(argv[i],"-s")==0)
                strcpy(data_file, argv[i+1]);
            else if(strcmp(argv[i],"-o")==0)
                strcpy(out_file, argv[i+1]);
            else if(strcmp(argv[i],"-c")==0)
                strcpy(solution_file, argv[i+1]);
            else if(strcmp(argv[i],"-t")==0)
                MAX_TIME = atoi(argv[i+1]);
        }
        //printf("data_file= %s, output_file= %s, sln_file=%s, max_time=%d", data_file, out_file, solution_file, MAX_TIME);
    }
    struct problem_struct** my_problems = load_problems(data_file);
    if(strlen(solution_file)<=0)
    {
        if(strcmp(out_file,"")==0) strcpy(out_file, "my_solutions.txt"); //default output
        FILE* pfile = fopen(out_file, "w"); //open a new file
        fprintf(pfile, "%d\n", num_of_problems); fclose(pfile);
        for(int k=0; k<num_of_problems; k++) 
        {
            best_sln.objective=0; best_sln.feasibility=0;
            for(int run=0; run<NUM_OF_RUNS; run++)
            {
                srand(RAND_SEED[run]);
                particleSwarm(my_problems[k]); //call PSO algorithm
            }
            output_solution(&best_sln,out_file); //output the solution
        }
    }
    for(int k=0; k<num_of_problems; k++)
    {
       free_problem(my_problems[k]); //free problem data memory
    }
    free(my_problems); //free problems array
    if(best_sln.x!=NULL && best_sln.cap_left!=NULL){ free(best_sln.cap_left); free(best_sln.x);} //free global
    return 0;
}