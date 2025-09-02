#Jarow Myers
import math
import random
import JMdatasets

print("                                     ")
print("-------------------------------------")
print("Welcome to Jarow's Genetic Algorithm for the P-Median Problem.")
print("-------------------------------------")
print("                                     ")

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

P = 10
if (problem_index == 0):
    print(f"Choose a Median-Set Size (for a graph of 10 points)")
elif(problem_index == 1):
    print(f"Choose a Median-Set Size (for a graph of 50 points)")
else:
    print(f"Choose a Median-Set Size (for a graph of 100 points)")
P = int(input())


print("Choose a selection technique:\n[1]Roulette     [2]Tournament")
selection_index = int(input())-1

print("Choose a crossover operator:\n[1]Double Point     [2]Uniform")
crossover_index = int(input())-1

print("Choose a crossover rate [0-1]")
CROSSOVER_RATE = float(input())

print("Choose a mutation operator:\n[1]Single Point     [2]Nearest Neighbor")
mutation_index = int(input())-1

print("Choose a mutation rate [0-1]")
MUTATION_RATE = float(input())


selection_technique = ["roulette", "tournament"]

#for fix-ups on crossover, prefer the new alleles to the old ones
crossover_operator = ["double_point", "uniform"]

#single_swap removes a random median and adds another
#nearest_neighbor takes a median and makes its closest neighboring point a median instead
mutation_operator = ["single_point", "nearest_neighbor"]

POPULATION_SIZE = 100

MUTATION = mutation_operator[mutation_index]

SELECTION = selection_technique[selection_index]

CROSSOVER = crossover_operator[crossover_index]

TOURNAMENT_RATE = .75
# if (SELECTION == "tournament"):
#     #get tournament rate
#     pass

print("\n\nStructure for GA:\n------------------------------\n")
print(f"Problem: {problem_type}")
print(f"Selection: {SELECTION}\nCrossover: {CROSSOVER}\nMutation: {MUTATION}")
print(f"Crossover Rate: {CROSSOVER_RATE}\nMutation Rate: {MUTATION_RATE}")
if (SELECTION == "tournament"):
    print(f"Tournament Rate: {TOURNAMENT_RATE}")
print(f"Population Size: {POPULATION_SIZE}\nGraph Size: {len(GRAPH)}\nMedian-Set Size: {P}")
print("\n------------------------------")


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

#returns a chromosome with a specified length and an empty median set
def createBlankChromosome(length = CHROMOSOME_SIZE):
    empty_list = []
    for i in range(length):
        empty_list.append(0)

    return Chromosome(empty_list)

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

#creates two children by crossing over two parents using the given technique
def crossover(parent1, parent2, technique = CROSSOVER):
    
    if (technique == "double_point"):
        point1 = random.randint(0,len(parent1.genes))
        point2 = random.randint(point1, len(parent1.genes))
        p1_genes = parent1.genes
        p2_genes = parent2.genes
        child1 = Chromosome(p1_genes[:point1] + p2_genes[point1:point2] + p1_genes[point2:])
        child2 = Chromosome(p2_genes[:point1] + p1_genes[point1:point2] + p2_genes[point2:])

        #fix child genes to have right number of medians
        fixGenes(child1, parent2)
        fixGenes(child2, parent1)

        return child1, child2

    elif(technique == "uniform"):
        #create uniform string
        uniform_string = []
        for i in range(len(parent1.genes)):
            if (random.random() < .5):
                uniform_string.append(0)
            else:
                uniform_string.append(1)
        child1_genes = []
        child2_genes = []
        for i in range(len(parent1.genes)):
            if uniform_string[i] == 0:
                child1_genes.append(parent1.genes[i])
                child2_genes.append(parent2.genes[i])
            else:
                child1_genes.append(parent2.genes[i])
                child2_genes.append(parent1.genes[i])
        child1 = Chromosome(child1_genes)
        child2 = Chromosome(child2_genes)

        #fix genes
        fixGenes(child1, parent2)
        fixGenes(child2, parent1)

        return child1, child2


#make new_child have the correct number of medians, preferring to use new_parent's medians
def fixGenes(new_child, new_parent, num_medians = P):
    child_medians = new_child.getMedians()
    parent_medians = new_parent.getMedians()
    
    #if child needs more medians
    if (len(child_medians) < num_medians):
        num_needed = num_medians - len(child_medians)
        
        #get required number of medians from parent
        potential_meds = set(parent_medians) - set(child_medians)
        new_meds = random.sample(potential_meds, num_needed)
        
        #add the medians to child
        for med_idx in new_meds:
            new_child.genes[med_idx] = 1
    
    #elif child has too many medians, randomly remove some
    elif (len(child_medians) > num_medians):
        #get required medians
        num_needed = len(child_medians) - num_medians
        to_delete = random.sample(child_medians, num_needed)

        #delete the medians from child
        for med_idx in to_delete:
            new_child.genes[med_idx] = 0




#alters the genes of an individual using the given mutation technique
def mutate(individual, technique=MUTATION):
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


    return


#creates a random population of chromosomes
def createInitialPopulation(pop_size = POPULATION_SIZE):
    init_pop = []
    for i in range(pop_size):
        init_pop.append(createRandomChromosome())
    return init_pop

#builds a parent_pool from a population using the given selection_technique
def buildParentPool(population, selection_technique = SELECTION):
    parent_pool = []
    if (selection_technique == "roulette"):
        #get adjusted fitnesses since fitness is to be minimized
        adjusted_fitnesses, adjusted_total = getAdjustedFitness(population)
        for i in range(len(population)):
            rdm_value = adjusted_total * random.random()
            total = 0
            for (index, value) in enumerate(adjusted_fitnesses):
                total += value
                if (total > rdm_value):
                    parent_pool.append(population[index])
                    break

    elif(selection_technique == "tournament"):

        for i in range(len(population)):
            contenders = random.sample(population, 2)
            if (fitness(contenders[0]) > fitness(contenders[1])):
                #add better individual TOURNAMENT_RATE% of time
                if (random.random() < TOURNAMENT_RATE):
                    parent_pool.append(contenders[0])
                else: #else add the worse individual
                    parent_pool.append(contenders[1])
            else:#else second is better
                #add better individual TOURNAMENT_RATE% of time
                if (random.random() < TOURNAMENT_RATE):
                    parent_pool.append(contenders[1])
                else: #else add the worse individual
                    parent_pool.append(contenders[0])

    return parent_pool

#creates a child pool from a given parent pool using the provided crossover technique and rate
def buildChildPool(parent_pool, crossover_technique = CROSSOVER, crossover_rate = CROSSOVER_RATE):
    children = []

    for i in range(len(parent_pool)//2):
        parent1, parent2 = random.sample(parent_pool, 2)
        if (random.random() < crossover_rate):
            child1, child2 = crossover(parent1, parent2, technique=crossover_technique)
        else:
            child1 = parent1
            child2 = parent2
        children.append(child1)
        children.append(child2)

    return children

#loops through a population and mutates individual based on mutation rate and technique
def mutatePopulation(population, mutation_rate = MUTATION_RATE, technique=MUTATION):
    for ind in population:
        if (random.random() < mutation_rate):
            mutate(ind, technique=technique)


#gets the top two individuals from a population
def getElites(population):
    elites = []

    #get first two, order by fitness
    #so first is better than second
    first = population[0]
    first_fitness = fitness(first)
    second = population[1]
    second_fitness = fitness(second)
    if(second_fitness<first_fitness):
        tmp = first
        tmpf = first_fitness
        first = second
        first_fitness = second_fitness
        second = tmp
        second_fitness = tmpf
    
    #loop through rest
    for ind in population[2:]:
        ind_fitness = fitness(ind)
        #if better than first, make it first and shift old first to second
        if(ind_fitness < first_fitness):
            second = first
            second_fitness = first_fitness
            first = ind
            first_fitness = ind_fitness

        #else if better than second, replace second
        elif(ind_fitness < second_fitness):
           second = ind
           second_fitness = ind_fitness

    elites.append(first)
    elites.append(second)
    return elites

#finds the individual with the best fitness in the given population
#RETURNS: the chromosome with the best (lowest) fitness
def findBestIndividual(population):
    best = population[0]
    best_fitness = fitness(best)
    for ind in population:
        if (fitness(ind) < best_fitness):
            best = ind
            best_fitness = fitness(ind)

    return best

#returns total fitness of a population
def getTotalFitness(population):
    total = 0
    for ind in population:
        total += fitness(ind)
    return total

#returns a list of adjusted fitnesses corresponding to each individual, and the sum
#adjusted fitness is so that a lower fitness value is valued higher
def getAdjustedFitness(population):
    adjusted_list = []
    total_fitness = getTotalFitness(population)
    for ind in population:
        new_fitness = total_fitness/fitness(ind)
        adjusted_list.append(new_fitness)
    adjusted_total = sum(adjusted_list)
    return adjusted_list, adjusted_total


############ Main Execution #############
working_population = createInitialPopulation()
num_iterations_no_change = 0
generations_run = 0
overall_best = findBestIndividual(working_population)
overall_best_fitness = fitness(overall_best)
print(f"Best individual at start:\n{overall_best.genes}")
print(f"Fitness: {overall_best_fitness}\n")
print("------------------------------")
while((generations_run < 1000) and (num_iterations_no_change < 100)):

    #get top two parents
    elites = getElites(working_population)

    #get parent pool using selection
    parent_pool = buildParentPool(working_population)

    #reproduce parents to get children
    child_population = buildChildPool(parent_pool)

    #mutate the children
    mutatePopulation(child_population)

    #put elites in final population
    child_population[0] = elites[0]
    child_population[1] = elites[1]

    #replace old population with new population
    working_population = child_population

    #find the most fit of the new population
    best_individual = findBestIndividual(working_population)
    best_individual_fitness = fitness(best_individual)

    #if the best of the new pop is better than any before, take note
    #if it is not better, take note of how long the old best has stood
    if (best_individual_fitness < overall_best_fitness):
        num_iterations_no_change = 0
        overall_best = best_individual
        overall_best_fitness = best_individual_fitness
    else:
        num_iterations_no_change += 1
    
    generations_run += 1
    #print an update every 20 generations
    if ((generations_run % 20) == 0 ):
        print(f"{generations_run} generations run.")
        print(f"Best individual so far:\n{overall_best.genes}")
        print(f"Fitness: {overall_best_fitness}")
        print("------------------------------")
#end while

print(f"Algorithm has finished after {generations_run} generations.")
print(f"Best individual found:\n{overall_best.genes}")
median_list = overall_best.getMedians()
results = []
for med in median_list:
    results.append(GRAPH[med])
print(f"Medians: {results}")
print(f"Fitness: {overall_best_fitness}")
