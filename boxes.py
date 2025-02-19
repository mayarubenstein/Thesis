import numpy
from sortedcontainers import SortedList

class radiusBox:

    # initialize the box with its boxCoords, an ennumeration of
    # the boxes in the volumeHT
    def __init__(self, number):
        self.boxNumber = number
        self.contents = SortedList()

    # returns the number of particles in the box
    def population(self):
        return len(self.contents)

    # add an item to the box
    def add(self, item):
        self.contents.add(item)

    # delete a specific item from the box
    def delete(self, item):
        self.contents.remove(item)

    # evaluates if the box is empty
    def isEmpty(self):
        return len(self.contents) == 0
    
class overlapBox:
    
    #initialize the box with its boxCoords, an ennumeration of the boxes
    #in volumeHT
    def __init__(self, number):
        self.boxNumber = number
        self.contents = {}
        self.size = 0
    
    #returns the number of particles in the box
    def population(self):
        return self.size
    
    #add an item to the box
    def add(self, item):
        self.contents[item] = item
    
    #delete an item from the box
    def delete(self, item):
        try:
            self.contents.pop(item)
        except:
            return


#like radius box class, but counts how many particles of each type
#assumes molecule classes species attribute
#also includes functionality for checking collisions
class augRadiusBox:
    # initialize the box with its boxCoords, an ennumeration of
    # the boxes in the volumeHT
    def __init__(self, number, listOfParticleTypes):
        self.boxNumber = number
        self.contents = SortedList()
        speciesCounter = {}
        for species in listOfParticleTypes:
            speciesCounter[species] = 0
        self.speciesCounter = speciesCounter

    # returns the number of particles in the box
    def population(self):
        return len(self.contents)

    # add an item to the box
    def add(self, item):
        self.contents.add(item)
        self.speciesCounter[item.species] +=1

    # delete a specific item from the box
    def delete(self, item):
        self.contents.remove(item)
        self.speciesCounter[item.species] -= 1

    # evaluates if the box is empty
    def isEmpty(self):
        return len(self.contents) == 0
    
    # return the population of a given species
    def speciesPopulation(self, species):
        try:
            return self.speciesCounter[species]
        except:
            raise Exception("No particles of species " + str(species) +
                            " in this simulation.")
    
    # is the move allowed?
    def allowedMove(self, item, coord):
        import numpy
        for particle in self.contents:
            if particle == item:
                continue
            elif numpy.norm(coord - particle.pos)<= item.rad + particle.rad:
                return False
        return True
    

#like overlap box, except keeps count of each species, assuming species 
# attribute in the particles
class augOverlapBox:
    
    #initialize the box with its boxCoords, an ennumeration of the boxes
    #in volumeHT
    def __init__(self, number, listOfParticleTypes):
        self.boxNumber = number
        self.contents = {}
        self.size = 0
        speciesCounter={}
        for species in listOfParticleTypes:
            speciesCounter[species] = 0
        self.speciesCounter = speciesCounter
    
    #returns the number of particles in the box
    def population(self):
        return self.size
    
    #add an item to the box
    def add(self, item):
        self.contents[item] = item
        self.size +=1
        self.speciesCounter[item.species] +=1
    
    #delete an item from the box
    def delete(self, item):
        try:
            self.contents.pop(item)
            self.size -=1
            self.speciesCounter[item.species] -=1
        except:
            return
    
    def speciesPopulation(self, species):
        try:
            return self.speciesCounter[species]
        except:
            raise Exception("No particles of species " + str(species) +
                            " in this simulation.")
