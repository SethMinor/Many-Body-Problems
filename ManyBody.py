# -*- coding: utf-8 -*-

"""
Simulating gravitational attraction between massive bodies, son!
The development of this horrific code can be blamed upon Seth Minor.

Compiling the code *should* produce an example of bounded motion in
an inverse-square potential.
"""

# imports
import pygame
import numpy as np
import sys

# initializations
pygame.init()

# parameters
width = 1400
height = 900
background_color = (31,34,38)
white = (255,255,255)
body_color = (188,196,209)
fps = 60
t = 0
dt = 1#1/fps
G = 0.000001   # universal graviational constant

# text stuff...
time_font = pygame.font.SysFont('euphemia',12)
info_font = pygame.font.SysFont('euphemia',12)

# the background window
window = pygame.display.set_mode((width,height))
pygame.display.set_caption("graverino")
window.fill(background_color)


# classes
class Body:
    # constructor and attributes
    def __init__(self, name, mass, color, x0, y0, vx_0, vy_0):
        self.name = name
        self.mass = mass
        self.color = color   # this should be in (r,g,b) format
        self.pos = np.array([x0, y0])
        self.vel = np.array([vx_0, vy_0])
        self.pos_tuple = (x0,y0)
        self.vel_tuple = (vx_0,vy_0)
    # get's
    def get_body_mass(self):
        return self.mass
    def get_body_color(self):
        return self.color
    def get_body_pos(self):
        return self.pos
    def get_body_vel(self):
        return self.vel
    def get_body_pos_tuple(self):
        return self.pos_tuple
    def get_body_vel_tuple(self):
        return self.vel_tuple


# Surfaces n' objects n' stuff
b_1 = Body("Body 1", 14.6, body_color, 180,100, -3,3)
b_2 = Body("Body 2", 15.2, body_color, 249,350, -1,2)
b_3 = Body("Body 3", 19.7, body_color, 607,756, -2,-1)
b_4 = Body("Body 4", 7.8, body_color, 0,14, 6,5)
b_5 = Body("Body 5", 28.3, body_color, 600,34, 3,1)
body_list = [b_1,b_2,b_3,b_4,b_5]


# Properties of the bodies
def radius(body):   # make radius proportional to mass
    k = 0.62        # scaling constant
    r = np.ceil(k * body.mass)
    return r

def print_body_stats(body):
    print("\n")
    print(body.name)
    print("Mass:",body.mass,"kg")
    print("Color:",body.color)
    print("Initial position:",body.pos)
    print("Initial velocity:",body.vel)
    
def print_body_list_stats(body_list):
    for body in body_list:
        print_body_stats(body)
    print("\n")


# PHYSICS functions, hamies
def get_distance(p1,p2):   # where p1=[x1,y1] and p2=[x2,y2]
    d = np.sqrt(np.power((p2[0]-p1[0]),2)+np.power((p2[1]-p1[1]),2))
    return d

def get_vector(p1,p2):   # returns a vector between points p1 and p2
    v = np.zeros((1,2))
    v[0,0] = p2[0]-p1[0]
    v[0,1] = p2[1]-p1[1]
    vprime = np.array([v[0,0],v[0,1]])
    return vprime

def get_unit_vector(p1,p2):   # returns the zero vector if p1==p2
    u = get_vector(p1,p2)
    d = get_distance(p1,p2)
    if (d != 0):              # don't divide by 0
        u = u/np.linalg.norm(u)
        return u
    else:
        return u

def get_forces_matrix(body_list):
    l = len(body_list)             # returns a matrix of forces
    F = np.zeros((l,l))            # where F_ij = |grav force pulling body_i towards body_j|
    i = 0; j = 0
    for body1 in body_list:
        for body2 in body_list:
            #print("for body1 = b_",i+1)
            #print("for body2 = b_",j+1)
            x1 = body1.pos[0]
            y1 = body1.pos[1] 
            x2 = body2.pos[0]
            y2 = body2.pos[1]
            if body1 == body2:
                F[i,j] = 0
                #print("force is zero")
                #print("F = ",F)
            elif body1 != body2:
                p1 = np.array([x1,y1])
                p2 = np.array([x2,y2])
                r = get_distance(p1,p2)
                F[i,j] = G *(body1.mass*body2.mass)/(np.power(r,2))   # governing force law, Newton's Law of Grav... son!
                #print("force is nonzero")
                #print("F = ",F)
            if (j+1) < len(body_list):
                #print("j = ",j,"  -->  j = ",j+1,"\n")
                j += 1
        if (i+1) < len(body_list):
            #print("i = ",i,"  -->  i = ",i+1,"\n")
            j = 0
            i += 1
    return F

def get_accelerations_matrix(body_list):   # returns a matrix of accelerations
    F = get_forces_matrix(body_list)       # where A_ij = |acceleration of body_i towards body_j|
    l = len(body_list)
    A = np.zeros((l,l))
    for i in range(l):
        for j in range(l):
            A[i,j] = F[i,j] / body_list[i].mass
    return A

# This function returns the cumulative acceleration vector on body_i at a certain time,
# where body_i is the i-th body in 'body_list'...
def get_acceleration_vector(i,body_list):          # returns a vector a_i = sum_over_j(a_ij), where |a_ij| = A_ij
    A = get_accelerations_matrix(body_list)        # and a_ij having base (xi,yi) and tip pointing towards (xj,yj)
    l = len(body_list)                             # (but really it's at the origin because vectors)
    a_i = np.array([0,0])
    p1 = np.array([body_list[i].pos[0],body_list[i].pos[1]])   # p1 = position of body_i
    #print("\np1 =",p1)
    #print("For A =\n",A)
    #print("and l =",l,":")
    for j in range(l):
        p2 = np.array([body_list[j].pos[0],body_list[j].pos[1]])   # p2 = position of body_j
        #print("   p2 =",p2)
        a_ij = A[i][j] * get_unit_vector(p1,p2)
        #print("   a_ij =",a_ij)
        a_i = np.add(a_i,a_ij)
        #print("     a_i ----> a_i =",a_i,"\n")
    return a_i

def do_physics(body_list, dt):
    i = 0
    for body in body_list:
        #print("For body = b_",i+1,"during a timestep of",dt,"sec:")
        
        # Update positions based on previos timestep velocities,
        # with approximation x_new = x_old + v*dt ...
        #print("positon = (",body.pos[0],",",body.pos[1],")")
        body.pos[0] += body.vel[0]#*dt
        body.pos[1] += body.vel[1]#*dt
        #print("      --> (",body.pos[0],",",body.pos[1],")")
        
        index = body_list.index(body)
        a_i = get_acceleration_vector(index,body_list)
        #print("acceleration vector =",a_i)
        
        # Update velocity vectors,
        # with approximation v_new = v_old + a*dt ...
        #print("velocity = (",body.vel[0],",",body.vel[1],")")
        body.vel[0] += np.ceil(a_i[0])#*dt
        body.vel[1] += np.ceil(a_i[1])#*dt
        #print("       --> (",body.vel[0],",",body.vel[1],")\n")
        
        i += 1


# Collision-checking func's
def handle_wall_collisions(body_list):   # based upon conservation of linear momentum, non-elastic impacts
    for body in body_list:
        if (body.pos[0]+radius(body) >= width) or (body.pos[0]-radius(body) <= 0):   # checking for left and right wall collisions
            body.vel[0] = -body.vel[0]
        if (body.pos[1]-radius(body) <= 0) or (body.pos[1]+radius(body) >= height):   # upper and lower wall collisions
            body.vel[1] = -body.vel[1]
        # handling the weird wall-jiggling effect...
        if (body.pos[0] >= width):   # right
            body.pos[0] = width-radius(body)
        if (body.pos[0] <= 0):   # left
            body.pos[0] = radius(body)
        if (body.pos[1] <= 0):   # up
            body.pos[1] = radius(body)
        if (body.pos[1] >= height):   # down
            body.pos[1] = height-radius(body)
            

def handle_body_collisions(body_list):   # also based on conservation of linear momentum
    for body_i in body_list:
        p1 = np.array([body_i.pos[0],body_i.pos[1]])
        for body_j in body_list:
            p2 = np.array([body_j.pos[0],body_j.pos[1]])
            d = get_distance(p1,p2)
            u_ij = get_unit_vector(p1,p2)
            u_ji = -u_ij
            if (d <= radius(body_i)+radius(body_j)):   # check for particle collisions
                #body_i.vel = pygame.math.Vector2.reflect(body_i.vel,u_ij)   # reflect about u vector
                #body_j.vel = pygame.math.Vector2.reflect(body_j.vel,u_ji)
                pass

def handle_collisions(body_list):   # body and wall collisions
    handle_wall_collisions(body_list)
    handle_body_collisions(body_list)


# Drawing functions
def draw_body(body_list):
    for body in body_list:
        r = radius(body)
        position_tuple = (body.pos[0],body.pos[1])
        pygame.draw.circle(window, body_color, position_tuple, r)
        
def draw_time(t):
    timer_no_jiggle = time_font.render("Time elapsed: ",1,white)
    timer = time_font.render(str(t) + "s",1,white)
    window.blit(timer_no_jiggle,(width-timer_no_jiggle.get_width()-45,\
                                  height-timer.get_height()))
    window.blit(timer,(width-timer.get_width(), height-timer.get_height()))
    
def draw_body_stats(body_list):   # draw info about bodies to the screen, only works for 3 bodies
    pos_info = info_font.render("Positions,   b1: "+str(b_1.pos)+"   b2: "+str(b_2.pos)\
                                +"   b3: "+str(b_3.pos),1,white)
    window.blit(pos_info, (10,height-85))
    vel_info = info_font.render("Velocities,   b1: "+str(b_1.vel)+ "   b2: "+str(b_2.vel)\
                                + "   b3: "+str(b_3.vel),1,white)
    window.blit(vel_info, (10,height-55))
    a1 = get_acceleration_vector(0,body_list)
    a2 = get_acceleration_vector(1,body_list)
    a3 = get_acceleration_vector(2,body_list)
    acc_info = info_font.render("Accelerations,   b1: "+str(a1)+"   b2: "+str(a2)\
                                +"   b3: "+str(a3),1,white)
    window.blit(acc_info, (10,height-25))
    
def draw_window(t,body_list):
    window.fill(background_color)
    draw_time(t)
    #draw_body_stats(body_list)
    draw_body(body_list)
    pygame.display.update()


# Exceptions
# def exceptions(body_list):
#     run = True
#     for body in body_list:
#         if (body.x < 0 or body.y < 0):
#             sys.exit("Exception: positions must be nonnegative values...")
#     return run


# Main animation loop
def main():
    clock = pygame.time.Clock()
    run = True
    
    print_body_list_stats(body_list)   
    
    while run:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False
                pygame.quit()
                sys.exit("End of simulation...")
                
        t = pygame.time.get_ticks() / 1000   # time (in seconds)
        clock.tick(fps)
        
        # exceptions(body_list)
        
        # do the physics, check collisions, then draw the window...
        do_physics(body_list,dt)
        handle_collisions(body_list)
        draw_window(t,body_list)
    main()


# main driver code
if __name__ == "__main__":
    main()