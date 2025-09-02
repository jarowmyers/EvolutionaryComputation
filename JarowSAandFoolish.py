#Jarow Myers
import math
import random
import JMdatasets

print("                                     ")
print("-------------------------------------")
print("Welcome to Jarow's Simulating Annealing and Hill-Climbing Algorithm for the P-Median Problem.")
print("-------------------------------------")
print("                                     ")

print("Which algorithm would you like to use?\n[1] Simulated Annealing  [2] Foolish Hill-Climbing")
FOOLISH = int(input()) - 1

print("Which data set would you like to use?\n[1]Toy Problem (10 points)    [2]Medium Problem (50 points)  [3]Large Problem 1 (100 points) [4]Large Problem 2 (100 points) [5]Large Problem 3 (100 clustered points)")
problem_index = int(input())-1

problems = [JMdatasets.toy_data, JMdatasets.med_data, JMdatasets.large_data1, JMdatasets.large_data2, JMdatasets.large_data3]
if (problem_index == 0):
    problem_type = "Toy Problem"
elif(problem_index == 1):
    problem_type = "Medium Problem"
elif(problem_index == 2):
    problem_type = "Large Problem 1"
elif(problem_index == 3):
    problem_type = "Large Problem 2"
else:
    problem_type = "Large Problem 3"

#Graph represenation
dataset = problems[problem_index]
GRAPH = []
GRAPH_SIZE = CHROMOSOME_SIZE = len(dataset)
GRAPH = dataset

print("Choose a mutation operator:\n[1]Single Point     [2]Nearest Neighbor")
mutation_index = int(input())-1

#single_swap removes a random median and adds another
#nearest_neighbor takes a median and makes its closest neighboring point a median instead
mutation_operator = ["single_point", "nearest_neighbor"]

PERTURBATION = mutation_operator[mutation_index]

P = 10
if (problem_index == 0):
    print(f"Choose a Median-Set Size (for a graph of 10 points)")
elif(problem_index == 1):
    print(f"Choose a Median-Set Size (for a graph of 50 points)")
else:
    print(f"Choose a Median-Set Size (for a graph of 100 points)")
P = int(input())


#CHROMOSOME DESIGN: BITSTRING
#With one allele for each node, 0 if in demand set, 1 if in median set
#length of chromosome is equal to |V|
#number of 1s is equal to P
#either make all chromosomes feasible (only allow P ones) or deal with infeasibles
#for each chromosome, keep track of indices of median set so we dont have to loop every time to find it
class Chromosome:
    
    def __init__(self, genes):
        self.genes = genes #allele for each node, either 0 or 1
        #self.medians = self.getMedians() #list of indices of the P nodes in the median set

    def getMedians(self):#returns a list of indices of median nodes
        med_list = []
        for allele in range(len(self.genes)):
            if self.genes[allele] == 1:
                med_list.append(allele)
        return med_list

#returns a random chromosome with a specified length and number of medians
def createRandomChromosome(length = CHROMOSOME_SIZE, num_medians = P):
    random_genes = []

    for i in range(length):
        random_genes.append(0)

    rand_meds = random.sample(range(length), num_medians)
    for i in rand_meds:
        random_genes[i] = 1

    return Chromosome(random_genes)


#returns fitness
def fitness(chromosome):
    #minimize the distance from median set to demand set
    #
    #FITNESS = SUM OVER EACH NODE: THE DISTANCE TO NEAREST MEDIAN NODE
    #
    fitness_value = 0

    for index, allele in enumerate(chromosome.genes):
        if allele == 0: #if node represented in that allele is in demand set
            fitness_value += distanceToMedian(index, chromosome)

    return fitness_value

#finds the distance between the node at allele_index and the nearest median node
def distanceToMedian(allele_index, chromosome):
    closest_median_dist = math.inf
    for med_index, median in enumerate(chromosome.getMedians()):#find the distance to closest median
        closest_median_dist = min(dist(GRAPH[allele_index], GRAPH[median]), closest_median_dist)
    return closest_median_dist

#Euclidean distance for 2d points
def dist(p, q):
    return math.sqrt(math.pow((p[0]-q[0]), 2) + math.pow((p[1]-q[1]), 2))


#alters the genes of an individual using the given mutation technique
def perturb(solution, technique=PERTURBATION):
    individual = Chromosome(solution.genes[0:])
    if (technique == "single_point"):
        #get medians
        ind_meds = individual.getMedians()
        #choose a random median
        rand_med_idx = ind_meds[random.randint(0, len(ind_meds) - 1)]
        
        #get a random non-median, set it to a median
        rand_idx_list = list(set(range(len(individual.genes))) - set(individual.getMedians()))
        rand_idx = rand_idx_list[random.randint(0, len(rand_idx_list) - 1)]
        individual.genes[rand_idx] = 1

        #make old median not a median
        individual.genes[rand_med_idx] = 0
    elif (technique == "nearest_neighbor"):
        #get medians
        ind_meds = individual.getMedians()
        #choose a random median
        rand_med_idx = ind_meds[random.randint(0, len(ind_meds) - 1)]
        #make it not a median
        individual.genes[rand_med_idx] = 0

        #find the nearest non-median point to the old median
        non_meds = list(set(range(len(individual.genes))) - set(individual.getMedians()))
        nearest_idx = non_meds[0]
        nearest_distance = dist(GRAPH[rand_med_idx], GRAPH[nearest_idx])
        for idx in non_meds:
            if(idx != rand_med_idx):
                point_dist = dist(GRAPH[idx], GRAPH[rand_med_idx])
                if (point_dist < nearest_distance):
                    nearest_idx = idx
                    nearest_distance = point_dist

        #make it a median
        individual.genes[nearest_idx] = 1


    return individual






#SA parameters
init_temp = 10
current_temp = init_temp
max_iterations = 1000
current_iterations = 0
temp_factor = .95
iteration_factor = 1.02

init_solution = createRandomChromosome(GRAPH_SIZE, P)
current_solution = Chromosome(init_solution.genes)

print("\n\nStructure for Algorithm:\n------------------------------\n")
print(f"Problem: {problem_type}")
print(f"Perturbation: {PERTURBATION}")
print(f"Graph Size: {len(GRAPH)}\nMedian-Set Size: {P}")
print("\n------------------------------")


total_iterations_run = 0
while (current_temp > 1):
    while(current_iterations < max_iterations):
        new_solution = perturb(current_solution)
        current_fitness = fitness(current_solution)
        new_fitness = fitness(new_solution)
        rdm = random.random()
        if (new_fitness < current_fitness):
            current_solution = new_solution
        elif (not FOOLISH):
            if (rdm < (math.pow(math.e, ((current_fitness - new_fitness)/current_temp)))):
                current_solution = new_solution
        
        # if(total_iterations_run % 100 == 0):
        #     print(current_fitness)

        current_iterations += 1
        total_iterations_run += 1
    current_temp = temp_factor * current_temp
    max_iterations = int(iteration_factor * max_iterations)
    current_iterations = 0
print("---------------------------------------")
print(f"Algorithm has finished after {total_iterations_run} iterations.")
print(f"Best individual found:\n{current_solution.genes}")
median_list = current_solution.getMedians()
results = []
for med in median_list:
    results.append(GRAPH[med])
print(f"Medians: {results}")
print(f"Fitness: {fitness(current_solution)}")
