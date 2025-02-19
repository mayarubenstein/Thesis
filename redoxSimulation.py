#constants
import numpy

ke = 1.44 #in eV * nm/e^2 #electric constant
r4 = 0.351 #nm #this is the radius of the -4 ions
r3 = 0.327 #nm #this is the radius of the -3 ions
kbt = 2568 #eV at 25C #Botzmann constant * temp
reducedDensity = 2e-3 * 1e-24 * 6.022e23 #density of solution in ions/nm**3
d3 = 7.2e-6 *1e14 #nm^2/sec #diffusion constant of -3 ion
d4 = 6.7e-6 *1e14 #nm^2/sec #diffusion constant of -4 ion
f = 1/kbt
alpha = 0.5
k = 1e-3 #s^-1 #this is the reaction coefficient
k = 1
E0 = 0.2 #eV/e #potential
timestep = 1e-11 #this is the timestep we are using here. For explanation why, see writeup

radius = {'3': r3, '4': r4} #radii for each ion
diffusionConstant = {'3': d3, '4': d4} #diffusion constant for each species

#molecule calss for type of sepcies
class molecule:
    def __init__(self, pos, species):
        self.pos = pos #the molecule's position
        self.species = species #the molecule's species; here given by charge
        self.rad = radius[species] #save the molecule's radius
        self.diffusion = diffusionConstant[species] #molecule's diffusion constant
        
    #function that simulates reducing the molecule from -3 to -4
    def reduce(self):
        if self.species != '3':  #check to make sure the particle has charge -3
            raise Exception("can't reduce")
        self.species = '4' #change to -4
        self.rad = radius['4'] #update radius
        self.diffusion = diffusionConstant[self.species] #update diffusion constant
    
    #simulate molecule being oxidized from -4 to -3
    def oxidize(self):
        if self.species != '4': #check the particle has charge -4
            raise Exception("can't oxidize")
        self.species = '3'  #change to -3
        self.rad = radius['3'] #update radius
        self.diffusion = diffusionConstant[self.species] #update diffusion constant
    
    #simulate reaction using Butler-Volmer formalism
    def react(self, potential):
        random = numpy.random.uniform()
        if self.species == "3":
            prob = k * numpy.exp((1- alpha)*f*(potential - E0)) #compute probability particle reacts based on rate constant
            if random < prob:
                self.reduce() #reduce if random number is less than prob
        elif self.species == "4": #same for oxidation
            prob = k * numpy.exp(-alpha *f * (potential - E0))
            if random < prob:
                self.oxidize()


def pbc(self, vector, lengths):
    newVector = numpy.zeros(len(vector))
    for i in range(len(vector)):
        coord = vector[i]
        while coord < -lengths[i]/2:
            newVector[i] = coord + lengths[i]
            coord += lengths[i]
        while coord > lengths[i]/2:
            newVector[i] = coord - lengths[i]
            coord -= lengths[i]
        else:
            newVector[i] = coord
    return newVector

class volumeHT:
    
    def __init__(self, x, y, z, c, d, boxType):
        self.conc = c #concentration
        self.lengths = tuple([x, y, z]) #save all dimensions for iterating
        tempS = numpy.cbrt(d/c) #compute the approximate box size using
        # the above formula (temp here is temporary not temperature)
        s = [] #save the actual box dimensions
        numberOfBoxes = []
        for length in self.lengths: #follow above formula to find
                                    # integer number of boxes in each dimension
            b = numpy.round(length/tempS)
            sDim = length/b #also find real box length
            numberOfBoxes.appen(b) #save these values
            s.append(sDim)
        self.numberOfBoxesPerDim = tuple(numberOfBoxes) #save to object variable
        self.sideLengths = tuple(s)
        self.boxDict = {}
        for i in range(self.numberOfBoxesPerDim[0]):
            for j in (self.numberOfBoxesPerDim[1]):
                for k in (self.numberOfBoxesPerDim[2]):
                    self.boxDict[(i, j, k)] = boxType((i, j, k))
    
    # find the box corresponding to the coordinates
    def coordsToBox(self, coords):
        boxCoords = []
        coords = coords + numpy.array([self.lengths[0]/2,
                                       self.lengths[1]/2, self.lengths[2]])
        for i in range(3):
            boxCoords.append(numpy.floor(coords[i]/self.sideLengths[i]))
        return numpy.array(boxCoords)
    
    #find the adjacent boxes
    def adjacent(self, box):
        lst = []
        for i in [-1, 0, 1]:
            new_x = box[0] + i
            if new_x < 0 or new_x >= self.numberOfBoxesPerDim[0]:
                continue
            for j in [-1, 0, 1]:
                new_y = box[1] + j
                if new_y < 0 or new_y >= self.numberOfBoxesPerDim[1]:
                    continue
                for k in [-1, 0, 1]:
                    new_z = box[2] + k
                    if new_z < 0 or new_z >= self.numberOfBoxesPerDim[2]:
                        continue
                    if (i, j, k) != (0, 0, 0):
                        lst.append((new_x, new_y, new_z))
        return lst

    # put a particle in the hash table
    def put(self, item, coords):
        box = self.coordsToBox(coords)
        self.boxDict[box].add(item)
    
    # deletes an item from the hash table
    def delete(self, item, coords):
        box = self.coordsToBox(coords)
        self.boxDict[box].delete(item)
    
    # get a list of all particles in the HashTable
    def contents(self):
        allItems = []
        allBoxes = self.boxDict.values()
        for box in allBoxes:
            for ob in box.contents:
                allItems.append(ob)
        return allItems