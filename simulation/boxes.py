import numpy
from sortedcontainers import SortedList

class radiusBox:

    # initialize the box with its boxCoords, an ennumeration of
    # the boxes in the volumeHT
    def __init__(self, number,  centerCoord, sideLengths, listOfParticleTypes = None):
        self.boxNumber = number
        self.contents = SortedList()
        self.center = centerCoord
        self.sideLengths = sideLengths

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
    
    # returns True if the particle is less than (s/2) - r form the
    # center and False otherwise, where s is the minimum side length
    # of the box
    def fromCenter(self, radius, coord):
        s = numpy.min(self.sideLengths)
        return numpy.norm(self.center - coord) < ((s/2) - radius)

    # returns a list of all the particles in the box
    def particles(self):
        return list(self.contents)
    
    # is the move allowed? Check other particles in box
    def allowedMoveInBox(self, item, coord):
        import numpy
        for particle in self.contents:
            if particle == item:
                continue
            elif numpy.norm(coord - particle.pos)<= item.rad + particle.rad:
                return False
        return True
    
    # is the move allowed? Check is particle from another box collides
    # with particles from this box
    def allowedMoveOutOfBox(self, coord):
        import numpy
        for particle in self.contents:
            if numpy.norm(coord - particle.pos)<= item.rad + particle.rad:
                return False
        return True
            
        
        
    
class overlapBox:
    
    #initialize the box with its boxCoords, an ennumeration of the boxes
    #in volumeHT
    def __init__(self, number, centerCoord = None, sideLengths = None, listOfParticleTypes = None):
        self.boxNumber = number
        self.contents = {}
        self.size = 0
        self.center = centerCoord
    
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
    
    # returns a list of all the particles in the box
    def particles(self):
        return self.contents.values()


#like radius box class, but counts how many particles of each type
#assumes molecule classes species attribute
#also includes functionality for checking collisions
class augRadiusBox:
    # initialize the box with its boxCoords, an ennumeration of
    # the boxes in the volumeHT
    def __init__(self, number, centerCoord, sideLengths, listOfParticleTypes = None,):
        self.boxNumber = number
        self.contents = SortedList()
        speciesCounter = {}
        for species in listOfParticleTypes:
            speciesCounter[species] = 0
        self.speciesCounter = speciesCounter
        self.center = centerCoord
        self.sideLengths = sideLengths

    # returns the number of particles in the box
    def population(self):
        return len(self.contents)
    
    # return the population of a given species
    def speciesPopulation(self, species):
        try:
            return self.speciesCounter[species]
        except:
            raise Exception("No particles of species " + str(species) +
                            " in this simulation.")

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
    
    def fromCenter(self, coord, radius):
        s = numpy.min(self.sideLengths)
        return numpy.norm(self.center - coord) < ((s/2) - radius)
    
    # is the move allowed?
    def allowedMoveInBox(self, item, coord):
        import numpy
        for particle in self.contents:
            if particle == item:
                continue
            elif numpy.norm(coord - particle.pos)<= item.rad + particle.rad:
                return False
        return True

    # returns a list of all the particles in the box
    def particles(self):
        return list(self.contents)
    
    
    # is the move allowed? Check other particles in box
    def allowedMoveInBox(self, item, coord):
        import numpy
        for particle in self.contents:
            if particle == item:
                continue
            elif numpy.norm(coord - particle.pos)<= item.rad + particle.rad:
                return False
        return True
    
    # is the move allowed? Check is particle from another box collides
    # with particles from this box
    def allowedMoveOutOfBox(self, coord):
        import numpy
        for particle in self.contents:
            if numpy.norm(coord - particle.pos)<= item.rad + particle.rad:
                return False
        return True
    

#like overlap box, except keeps count of each species, assuming species 
# attribute in the particles
class augOverlapBox:
    
    #initialize the box with its boxCoords, an ennumeration of the boxes
    #in volumeHT
    def __init__(self, number, centerCoord = None, sidelengths = None, listOfParticleTypes = None):
        self.boxNumber = number
        self.contents = {}
        self.size = 0
        self.center = centerCoord
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
    
    # returns a list of all the particles in the box
    def particles(self):
        return self.contents.values()