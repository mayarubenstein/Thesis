import numpy
from sortedcontainers import SortedList
import copy


class porousElectrode:
    
    X = 0
    Y = 1
    Z = 2
    ENTER = 0
    COLLIDE = 1
    NOPORE = 2
    
    def __init__(self, poreDict, depth, wallStick, floorStick, surfaceStick, wallReact, floorReact, surfaceReact, wallTheta, floorTheta, surfaceTheta, timestep, surfaceProx, floorProx, surfaceParameters):
        for theta in [floorTheta, wallTheta, surfaceTheta]:
            assert theta>=0 and theta < numpy.pi/2, "theta must be between 0 and pi/2"
        for stick in [floorStick, wallStick, surfaceStick]:
            assert stick>=0 and stick <= 1, "stick must be in [0, 1]"
        self.pores = []
        for center, radius in poreDict.items():
            self.pores.append(pore(center, radius, depth, wallStick, floorStick, floorProx, wallReact, floorReact, wallTheta, floorTheta, timestep))
        self.aht = areaHT(*surfaceParameters, self.pores)
        self.stuck = SortedList()
        self.surfaceStick = surfaceStick
        self.surfaceTheta = surfaceTheta
        self.surfaceReact = surfaceReact
        self.surfaceProx = surfaceProx
        self.surfaceStuck = []
        
    # returns 0 if enters pore, 1 if enters pore but collision, and 2 if
    # it doesn't enter pore
    def entersPore(self, particle, newPos):
        if newPos[self.Z] - particle.rad > 0:
            return self.NOPORE
        rect = self.aht.getBox(newPos)
        possiblePores = self.aht.getPores(rect)
        for p in possiblePores:
            if p.entryExitCondition(particle, newPos):
                if p.collisionChecker(particle, newPos):
                    p.place(particle, newPos)
                    return self.ENTER
                return self.COLLIDE
        
    def movingAndShakingPores(self, potential):
        exit = []
        for p in self.pores:
            exit += p.movingAndShakingFloor(potential)
            exit += p.movingAndShakingWalls(potential)
            exit += p.movingAndShakingInterior(potential)
        return exit
            
    def movingAndShakingSurface(self, potential):
        released = []
        for particle in copy.copy(self.surfaceStuck):
            if self.surfaceReact:
                particle.react(potential)
            if self.surfaceRelease(particle):
                released.append(particle)
        return released
        
    def hit(self, particle, potential):
        newPosZ = particle.pos[self.Z]
        distCutoff = numpy.max(particle.rad, self.surfaceProx)
        if self.surfaceReact:
            if newPosZ <= distCutoff:
                particle.react(potential)
        if newPosZ <= particle.rad:
            self.stuck.add(particle)
            return True
        return False
            
    def Bernoulli(p):
        return numpy.random.binomial(1, p)
    
    def surfaceRelease(self, particle):
        leaving = self.Bernoulli(1 - self.surfaceStick)
        if leaving:
            self.surfaceStuck.remove(particle)
            return True
        return False
    
    def rebound(self, particle):
        angle = numpy.uniform(low = 0, high = self.theta)
        phi = numpy.uniform(low = 0, high = 2* numpy.pi)
        sigma = numpy.sqrt(2*particle.diffusion*self.timestep)
        rho = numpy.random.normal(0, sigma)
        xMove = rho * numpy.sin(angle) * numpy.cos(phi)
        yMove = rho * numpy.sin(angle) * numpy.sin(phi)
        zMove = rho * numpy.sin(angle)
        displace = numpy.array([xMove, yMove, zMove])
        return particle.pos + displace
    
    def interception(self, particle, newPos):
        radius = particle.rad
        if (newPos[self.Z] - radius > 0): # no interception in this case
            return newPos
        oldPos = particle.pos
        ak = oldPos[2]
        bk = newPos[2]
        factor = (-ak/(bk - ak))
        displacementVector = newPos-oldPos
        return displacementVector * factor + oldPos - displacementVector* radius/numpy.linalg.norm(displacementVector)

class pore:
    
    X = 0
    Y = 1
    Z = 2
    
    def __init__(self, center, radius, depth, wallStick, floorStick, floorProx, wallReact, floorReact, wallTheta, floorTheta, timestep):
        self.center = center
        self.radius = radius
        self.wallStuck = SortedList()
        self.floorStuck = SortedList()
        self.inPore = SortedList()
        self.moveable= SortedList()
        self.depth = depth
        self.floorStick = floorStick
        self.wallStick = wallStick
        self.floorProx = floorProx
        self.wallReact = wallReact
        self.floorReact = floorReact
        assert floorTheta>=0 and floorTheta < numpy.pi/2, "theta must be between 0 and pi/2"
        assert wallTheta>=0 and wallTheta < numpy.pi/2, "theta must be between 0 and pi/2"
        self.wallTheta = wallTheta
        self.floorTheta = floorTheta
        self.timestep = timestep
    
    # assume particle.pos = oldPos
    def entryExitCondition(self, particle, newPos):
        xy = self.xyPoreAtZ(particle, newPos, 0)
        return self.inPore(particle, xy)
    
    # find intersection point between the line from particle.pos to newPos
    # and the plane z = zPos 
    def posAtZ(self, particle, newPos, zPos):
        oldPos = particle.oldPos
        factor = (zPos-oldPos[self.Z])/(newPos[self.Z] - oldPos[self.Z])
        return self.coordsFromT(oldPos, newPos, factor)
            
    def xyInPore(self, particle, coord):
        return (coord[self.X] - self.center[self.X])**2 + (coord[self.Y] - self.center[self.Y])**2 <= (self.radius - particle.rad)**2
    
    # return intercept of the particle's path and the cylinder's walls
    def xyIntercept(self, oldPos, newPos):
        ax = oldPos[0]
        ay = oldPos[1]
        bx = newPos[0]
        by = newPos[1]
        cx = self.center[0]
        cy = self.center[1]
        A = (bx - ax)**2 + (by - ay)**2
        B = 2*((bx - ax)*(ax - cx) + (by - ay)*(ay - cy))
        C = (ax - cx)**2 + (ay - cy)**2 - self.radius**2
        sqTerm = numpy.sqrt(B*B - 4*A*C)
        t = (-B + sqTerm)/(2*A)
        if t < 0 or t > 1:
            return False
        return t
    
    def coordsFromT(self, oldPos, newPos, t):
        return (newPos - oldPos) * t + oldPos
        
        
    # assume particle.pos = oldPos
    def movingAndPlacing(self, particle, newPos):
        # check if particle is in the interior of the pore, not stuck to anything
        # (assume we already know that the particle is in the pore)
        # this time we know particle was in the pore before
        if self.xyInPore(self, particle, newPos) and (newPos[self.Z] > -self.depth + particle.rad):
            if self.collisionChecker(particle, newPos):
                particle.pos = newPos
                if particle not in self.moveable:
                    self.moveable.add(particle)
                return True
        # check if it hits the floor before it hits the walls
        if self.floorCollision(particle, newPos):
            return True
        # if hits walls
        if self.wallCollision(particle, newPos):
            return True
        return False
    
    def entryPlacement(self, particle, newPos):
        # check if particle is in the interior of the pore, not stuck to anything
        # (assume we already know that the particle is in the pore)
        # this time we alrady know it passes the collision test
        if self.xyInPore(self, particle, newPos) and (newPos[self.Z] > -self.depth + particle.rad):
            particle.pos = newPos
            self.inPore.add(particle)
            if particle not in self.moveable:
                self.moveable.add(particle)
            return
        # check if it hits the floor before it hits the walls
        if self.floorCollision(particle, newPos):
            return
        # if hits walls
        if self.wallCollision(particle, newPos):
            return
        
    def floorCollision(self, particle, newPos):
        oldPos = particle.pos 
        floorPos = -self.depth
        interceptedPos = self.posAtZ(oldPos, newPos, floorPos + particle.rad)
        if self.xyInPore(self, particle, interceptedPos):
            if self.collisionChecker(particle, interceptedPos): # check if collision
                particle.pos = interceptedPos # update position
                if particle not in self.inPore: # add to inpore if not already in the list
                    self.inPore.add(particle)
                self.floorStuck.add(particle)
                self.moveable.discard(particle)
                return True
        return False # if it collides, or doesn't pass
        # interception test, return False
    
    def wallCollision(self, particle, newPos):
        oldPos = particle.poos
        t = self.xyIntercept(oldPos, newPos)
        interceptedPos = self.coordsFromT(oldPos, newPos, t)
        assert t>= 0 and t<= 1
        if interceptedPos[self.Z] < -self.depth:
            return False
        if self.collisionChecker(particle, interceptedPos):
            particle.pos = interceptedPos
            if particle not in self.inPore:
                self.inPore.add(particle)
            self.wallStuck.add(particle)
            self.moveable.discard(particle)
            return True
        return False
        
        
    def collisionChecker(self, particle, position):
        for secondParticle in self.inPore:
            if secondParticle == particle:
                continue
            elif numpy.linalg.norm(position - secondParticle.pos) < particle.rad + secondParticle.rad:
                return False
        return True
    
    def Bernoulli(p):
        return numpy.random.binomial(1, p)
    
    def wallRelease(self, particle):
        leaving  = self.Bernoulli(1 - self.wallStick)
        if leaving:
            self.wallStuck.remove(particle)
            return True
        return False
    
    def floorRelease(self, particle):
        leaving = self.Bernoulli(1 - self.floorStick)
        if leaving:
            self.floorStick.remove(particle)
            return True
        return False
    
    def movingAndShakingFloor(self, potential):
        exit = []
        for particle in copy.copy(list(self.floorStuck)):
            if self.floorReact:
                particle.react(potential)
            if self.floorRelease(particle):
                newPos = self.reboundFloor(particle)
                if newPos[self.Z] > 0 and self.entryExitCondition(particle, newPos): # check if exiting pore
                    exit.append(particle)
                    self.inPore.discard(particle)
                else:
                    self.place(particle, newPos)
        return exit
    
    def movingAndShakingWalls(self, potential):
        exit = []
        for particle in copy.copy(list(self.wallStuck)):
            if self.floorReact:
                particle.react(potential)
            if self.wallRelease(particle):
                newPos = self.reboundWall(particle)
                if newPos[self.Z] > 0 and self.entryExitCondition(particle, newPos): # check if exiting pore
                    exit.append(particle)
                    self.inPore.discard(particle)
                else:
                    self.place(particle, newPos)
        return exit
    
    def movingAndShakingInterior(self, potential):
        exit = []
        for particle in copy.copy(list(self.moveable)):
            # see if it is close enough to floor or walls to react
            if self.floorReact:
                if particle.pos[self.Z] + self.depth < self.floorProx:
                    particle.react(potential)
            sigma = numpy.sqrt(2*particle.diffusion*self.timestep)
            moveVector = numpy.random.normal(0, sigma, 3)
            newPos = particle.pos + moveVector
            if newPos[self.Z] > 0 and self.entryExitCondition(particle, newPos): # if it is leaving the pore
                exit.append(particle)
            else:
                self.place(particle)
        return exit # what to do about collision checking here??? maybe just 
        # leave as limitation
        
    
    def bounce(self, theta, particle):
        angle = numpy.uniform(low = 0, high = theta)
        phi = numpy.uniform(low = 0, high = 2* numpy.pi)
        sigma = numpy.sqrt(2*particle.diffusion*self.timestep)
        rho = numpy.random.normal(0, sigma)
        xMove = rho * numpy.sin(angle) * numpy.cos(phi)
        yMove = rho * numpy.sin(angle) * numpy.sin(phi)
        zMove = rho * numpy.sin(angle)
        return numpy.array([xMove, yMove, zMove])
    
    def reboundFloor(self, particle):
        displace = self.bounce(self.floorTheta, particle)
        return particle.pos + displace
    
    def reboundWall(self, particle):
        mover = self.bounce(self.wallTheta, particle)
        x = mover[0]
        y = mover[1]
        z = mover[2]
        px = particle.pos[0] - self.center[0]
        py = particle.pos[1] - self.center[1]
        r = self.radius - particle.rad
        first = x*py/r - z*px/r
        second = -x*px/r -z*py/r
        displace = numpy.array([first, second, y])
        return particle.pos + displace

        
        
class areaHT:
    
    X = 0
    Y = 1
    Z = 2
    
    def __init__(self, sidelengths, spaceCoordinates, numberOfBoxesPerDim, pores):
    
        # presume that pores is an actual list of pores
        self.poreFinder = {}
        self.sidelengths = sidelengths
        self.spaceCoordinates = [spaceCoordinates[self.X],
                                       self.spaceCoordinates[self.Y]]
        for i in range(numberOfBoxesPerDim[self.X]):
            for j in range(numberOfBoxesPerDim[self.Y]):
                self.poreFinder[(i, j)] = []
        for pore in pores:
            boxes = self.findAllBoxes(pore.center, pore.radius)
            for box in boxes:
                self.poreFinder[box].append(pore)
    
    def xToBox(self, coord):
        coord = coord - self.spaceCoordinates[self.X]
        return numpy.floor(coord/self.sidelengths[self.X])
    
    def yToBox(self, coord):
        coord = coord - self.spaceCoordinates[self.Y]
        return numpy.floor(coord/self.sidelengths[self.Y])
    
    def findAllBoxes(self, center, radius):
        boxNums = []
        xNumMin = self.xToBox(center[0] - radius)
        xNumMax = self.xToBox(center[0] + radius)
        yNumMin = self.yToBox(center[1] - radius)
        yNumMax = self.yToBox(center[1] + radius)
        for i in range(xNumMin, xNumMax + 1):
            for j in range(yNumMin, yNumMax + 1):
                boxNums.append((i, j))
    
    def getBox(self, position):
        xPos = position[self.X]
        yPos = position[self.Y]
        return (self.xToBox(xPos), self.yToBox(yPos))
    
    def getPores(self, box):
        return self.poreFinder[box]
        