#constants
import numpy
import boxes

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



class particle:
    
    def __init__(self, pos):
        self.pos = pos
        
    # comparator (by address)
    def __lt__(self, other):
        return id(self) < id(other)

    def __eq__(self, other):
        return id(self) == id(other)
    
        

#molecule calss for type of sepcies
class molecule:
    
    radius = {'3': r3, '4': r4} #radii for each ion
    diffusionConstant = {'3': d3, '4': d4} #diffusion constant for each species

    def __init__(self, pos, species, reactorTracker = None):
        self.pos = pos #the molecule's position
        self.species = species #the molecule's species; here given by charge
        self.rad = self.radius[species] #save the molecule's radius
        self.diffusion = self.diffusionConstant[species] #molecule's diffusion constant
        if reactorTracker is not None:
            self.hasReactorChecker = True
            self.reactorTracker = reactorTracker
        else:
            self.hasReactorChecker = False

    
    #function that simulates reducing the molecule from -3 to -4
    def reduce(self):
        if self.species != '3':  #check to make sure the particle has charge -3
            raise Exception("can't reduce")
        self.species = '4' #change to -4
        self.rad = self.radius['4'] #update radius
        self.diffusion = self.diffusionConstant[self.species] #update diffusion constant
    
    #simulate molecule being oxidized from -4 to -3
    def oxidize(self):
        if self.species != '4': #check the particle has charge -4
            raise Exception("can't oxidize")
        self.species = '3'  #change to -3
        self.rad = self.radius['3'] #update radius
        self.diffusion = self.diffusionConstant[self.species] #update diffusion constant
    
    #simulate reaction using Butler-Volmer formalism
    def react(self, potential):
        if self.species == "3" and potential >= E0:
            self.reduce() #reduce if random number is less than prob
            if self.hasReactorChecker:
                self.reactorTracker.reduce()
        elif self.species == "4" and potential <= E0: #same for oxidation
            self.oxidize()
            if self.hasReactorChecker:
                self.reactorTracker.oxidize()
    
    # comparator (by address)
    def __lt__(self, other):
        return id(self) < id(other)

    def __eq__(self, other):
        return id(self) == id(other)
    
    def __hash__(self):
        return hash(id(self))
    
    
class reactorTracker:
    
    def __init__(self):
        self.oxidation = 0
        self.reduction = 0
        self.oxidationLog = []
        self.reductionLog = []
        self.alphCurrent = []
        
    def reduce(self):
        self.reduction +=1
    
    def oxidize(self):
        self.oxidation +=1
    
    def timestep(self):
        self.oxidationLog.append(self.oxidation)
        self.reductionLog.append(self.reduction)
        self.alphCurrent.append(self.oxidation - self.reduction)
        self.oxidation = 0
        self.reduction = 0
    
    

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
    
    def __init__(self, spaceCoordinates, c, d, boxType, listOfParticleTypes = None):
        self.boxType = boxType
        self.conc = c #concentration
        # spaceCoordinates should give the coordinates of the overall
        # volume in space. Should be in the form [xmin, xmax, ymin, ymax,
        # zmin, zmax]
        self.spaceCoordinates = tuple(spaceCoordinates)
        lengths = []
        for i in range(3):
            lengths.append(spaceCoordinates[2*i + 1] - spaceCoordinates[2*i])
        self.lengths = tuple(lengths) #save all dimensions for iterating
        tempS = numpy.cbrt(d/c) #compute the approximate box size using
        # the above formula (temp here is temporary not temperature)
        s = [] #save the actual box dimensions
        numberOfBoxes = []
        for length in self.lengths: #follow above formula to find
                                    # integer number of boxes in each dimension
            b = int(numpy.round(length/tempS))
            if b == 0:
                b = 1
            sDim = length/b #also find real box length
            numberOfBoxes.append(b) #save these values
            s.append(sDim)
        self.numberOfBoxesPerDim = tuple(numberOfBoxes) #save to object variable
        self.sideLengths = tuple(s) #these are the box lengths
        self.density = self.conc * self.lengths[0]*self.lengths[1]*self.lengths[2]/(self.numberOfBoxesPerDim[0]*self.numberOfBoxesPerDim[1]*self.numberOfBoxesPerDim[2])
        self.boxDict = {}
        for i in range(self.numberOfBoxesPerDim[0]):
            xcenter = spaceCoordinates[0*2] + i*self.sideLengths[0] + self.sideLengths[0]/2
            for j in range(self.numberOfBoxesPerDim[1]):
                ycenter = spaceCoordinates[1*2] + j*self.sideLengths[1] + self.sideLengths[1]/2
                for k in range(self.numberOfBoxesPerDim[2]):
                    zcenter = spaceCoordinates[2*2] + k*self.sideLengths[2] + self.sideLengths[2]/2
                    self.boxDict[(i, j, k)] = boxType((i, j, k), (xcenter, ycenter, zcenter), self.sideLengths, listOfParticleTypes)
    
    # find the box corresponding to the coordinates
    def coordsToBox(self, coords):
        boxCoords = []
        # subtract the minimum (x, y, z) in the box, beaucse boxes start counting at zero
        coords = coords - numpy.array([self.spaceCoordinates[0],
                                       self.spaceCoordinates[2], self.spaceCoordinates[4]])
        for i in range(3):
            boxCoords.append(numpy.floor(coords[i]/self.sideLengths[i]))
        return self.boxDict[tuple(boxCoords)] #everything has to happen in tuple; numpy
    # arrays are not hashable
    
    # is a given coordinate in this box
    def inBox(self, coord):
        for i in range(3):
            if not (coord[i] >= self.spaceCoordinates[2*i] and coord[i] <= self.spaceCoordinates[2*i + 1]):
                return False
        return True
    
    def getSurfaceParameters(self):
        return [self.sideLengths, self.spaceCoordinates, self.numberOfBoxesPerDim]
    
    #find the adjacent boxes
    def adjacent(self, box):
        lst = []
        boxNum = box.boxNumber
        for i in [-1, 0, 1]:
            new_x = boxNum[0] + i
            if new_x < 0 or new_x >= self.numberOfBoxesPerDim[0]:
                continue
            for j in [-1, 0, 1]:
                new_y = boxNum[1] + j
                if new_y < 0 or new_y >= self.numberOfBoxesPerDim[1]:
                    continue
                for k in [-1, 0, 1]:
                    new_z = boxNum[2] + k
                    if new_z < 0 or new_z >= self.numberOfBoxesPerDim[2]:
                        continue
                    if (i, j, k) != (0, 0, 0):
                        lst.append(self.boxDict[(new_x, new_y, new_z)])
        return lst 
    
    # put a particle in the hash table, check for collisions
    def put(self, item):
        coord = item.pos
        if not self.inBox(coord):
            return False
        box = self.coordsToBox(coord)
        # we don't need to check collision if in bulk region
        if (self.boxType == boxes.overlapBox) or (self.boxType == boxes.augOverlapBox):
            box.add(item)
            return True
        # must check for collision
        # first check if overlaps with any other particles in the box
        allowed = box.allowedMoveInBox(item, coord)
        if not allowed:
            return False
        # if the particle is close enough to the center of the box,
        # particles in adjacent boxes need not be checked
        if box.fromCenter(coord, item.rad):
            box.add(item)
            return True
        # if the particle is outside the box's center orbit, we must
        # check the particles in adjacent boxes for collision
        adjacentBoxes = self.adjacent(box)
        for b in adjacentBoxes:
            if not b.allowedMoveOutOfBox(item, coord):
                return False
        # if it does not collide with the partilces in adjacent boxes,
        # move the particle and return True
        box.add(item)
        return True
        
    
    # deletes an item from the hash table
    def delete(self, item):
        box = self.coordsToBox(item.pos)
        box.delete(item)
    
    # get a list of all particles in the HashTable
    def contents(self):
        allItems = []
        allBoxes = self.boxDict.values()
        for box in allBoxes:
            for ob in box.contents:
                allItems.append(ob)
        return allItems
    
    # attempt to move the particle item to the coordinate newPos
    # particle is deleted and readded for rehashing
    # if succeeds return True; otherwise return False
    def attemptMove(self, item, newPos):
        # save the old coord in case the move fails (ie collision)
        if not self.inBox(item.pos):
            return False
        oldPos = item.pos
        # delete the particle for rehashing
        try:
            self.delete(item)
        except:
            pass
        # update the position
        item.pos = newPos
        # rehash, checking for collisions
        if self.put(item):
            return True
        # if fails, return to old coordinates and readd to box
        item.pos = oldPos
        oldBox = self.coordsToBox(oldPos)
        # don't bother with checking for collisions with whole put()
        # regimen if just undoing
        oldBox.add(item)
        return False