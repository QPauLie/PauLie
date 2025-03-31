import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation
from paulie.helpers.recording import RecordGraph
import numpy as np
import math



def plot_graph(vertices, edges, edge_labels = None):
   graph = nx.Graph()
   graph.add_nodes_from(vertices)
   graph.add_edges_from(edges)
   pos = nx.spring_layout(graph)
   if edge_labels is not None:
       nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=edge_labels,font_color="red")
   nx.draw_networkx(graph, pos=pos)
   plt.show()


def plot_graph_by_nodes(nodes, commutators=[]):
    vertices, edges, edge_labels = nodes.get_graph(nodes, commutators)
    return plot_graph(vertices, edges, edge_labels)



def animation_graph(record: RecordGraph, interval=1000, repeat=False, storage=None):
   graph = nx.Graph()
   fig, ax = plt.subplots(figsize=(6,4))
   def clear():
       ax.clear()
       graph.remove_nodes_from(list(n for n in graph.nodes)) 


   def build_positions(edges, center, lighting):
       legs = []
       positions = {}

       for edge in edges:
           if center in edge:
               v = edge[1] if center == edge[0] else edge[0]
               legs.append([v])
       for i, leg in enumerate(legs):
           current = leg[0]
           while True:
               is_found = False
               for edge in edges:
                   v = edge[1] if current == edge[0] else edge[0]
                   if current in edge and v not in leg and v != center:
                       leg.append(v)
                       current = v
                       is_found = True
                       break
               if not is_found:
                   break
       legs.sort(key = lambda x: len(x))
       len_line = len(legs)

       positons = {}
       max_line = 7
       y_dist = 0.25
       pos_y = 0
       y = pos_y
       center_x = 0
       x_position_lighting = 0
       x_first = 0
       x_last = 0

       if len(legs) > 1:
           dist = 2.0/(8 if len_line > 7 else len_line)
           len_line = len(legs[len(legs)-1]) + len(legs[len(legs)-2]) + 1
           if len_line > max_line:
               len_line = max_line
           n = 0
           x = 1 + dist/2

           for v in reversed(legs[len(legs)-2]):
              
               x -= dist
               x_first = x if x_first == 0 else x_first
               positions[v] = np.array([x, y])
               n += 1

           x -= dist
           positions[center] = np.array([x, y])
           center_x = x
           n += 1
           direction = -1
           for v in legs[len(legs)-1]:
               if n > max_line:
                   x_last = x if x_last == 0 else x_last
                   x += direction * dist
                   y -= y_dist 
                   n = 0
                   max_line = 5
                   direction *= -1
                   x_position_lighting = (x + 1 - dist/2)/2 if x_position_lighting == 0 else x_position_lighting

               x += direction * dist
               positions[v] = np.array([x, y])
               n += 1

           x_last = x if x_last == 0 else x_last
           x_position_lighting = (x_first + x_last)/2 if x_position_lighting == 0 else x_position_lighting
           del legs[len(legs)-1]
           del legs[len(legs)-1]
           # other legs
           direction = 1
           if len(legs) > 0:
               n = len(legs)
               ang = 3*math.pi/(2*n)
               if ang >= math.pi/2:
                   ang = ang / 2
               c_ang = ang
               
               for leg in legs:
                   v = leg[0]
                   y = pos_y + dist*math.sin(c_ang)
                   x = center_x + dist*math.cos(c_ang)
                   positions[v] = np.array([x, y])
                   if len(leg) > 1:
                       v = leg[1]
                       y = pos_y + 2*dist*math.sin(c_ang)
                       x = center_x + 2*dist*math.cos(c_ang)
                       positions[v] = np.array([x, y])

                   c_ang += direction*ang
                   if c_ang > 3*math.pi/4:
                       direction *= -1 
                       c_ang = direction*ang
       elif len(legs) == 1:
             positions[legs[0][0]] = np.array([0, y])
             positions[center] = np.array([0.25, y])
             x_position_lighting = 0.125

       else:
             positions[center] = np.array([0, y])
             x_position_lighting = 0


       return positions, x_position_lighting


   def update(num):
       clear()
       frame = record.get_frame(num)
       ax.set_title(f"{frame.get_title()}")
       vertices, edges, edge_labels = record.get_graph(num)
       center = None
       with_labels = True

       if len(vertices) > 0:
           center = vertices[0]
           if len(center) > 10:
               edge_labels = None
               with_labels = False


       lighting = frame.get_lighting()

       if lighting:
           if len(lighting) > 10:
               edge_labels = None
               with_labels = False
           vertices.append(lighting)
           if vertices is not None:
               for v in vertices:
                   if frame.get_is_lits(v):
                       edges.append((lighting, v))
           if frame.get_is_dependent(lighting):
               dependent = lighting
               lighting = f"dependend {lighting}"
               edges.append((dependent, lighting))

       if frame.is_appending() and frame.is_removing():
           for v in vertices:
               if frame.get_is_removing(v):
                   vertices.remove(v)


       graph.add_nodes_from(vertices)
       graph.add_edges_from(edges)
       if frame.get_init():
           positions = nx.spring_layout(graph)
       else:
           if not record.get_is_prev(num):
               positions, x_position_lighting = build_positions(edges, center, lighting)
               record.set_positions(positions)
               record.set_x_position_lighting(x_position_lighting)
               positions[lighting] = np.array([x_position_lighting, 1])
           else:
               positions = record.get_positions()
               x_position_lighting = record.get_x_position_lighting()
               positions[lighting] = np.array([x_position_lighting, 1])


       #positions = nx.spring_layout(graph)

       color_map = []
       for node in graph:
           if lighting == node:
               if frame.get_is_dependent(node):
                   color_map.append('#2F4F4F')
               else:
                   color_map.append('red')
           else:
               if not frame.get_is_lits(node):
                   if frame.get_is_q(node):
                       color_map.append('#FF00FF')
                   else:
                       if frame.get_is_removing(node):
                           color_map.append('black')
                       else:
                           color_map.append('#cccccc')

               else:
                   if frame.get_is_appending(node):
                       color_map.append('#00FF00')
                   elif frame.get_is_contracting(node):
                       color_map.append('#008080')
                   else:
                       if frame.get_is_p(node):
                           color_map.append('#6A5ACD')
                       else:
                           if frame.get_is_replacing(node):
                               color_map.append('#8B008B')
                           else:
                               color_map.append('cyan')

       plt.axis("off")

       if edge_labels is not None:
           nx.draw_networkx_edge_labels(graph, pos=positions, edge_labels=edge_labels,font_color='red', hide_ticks=True, node_size=60,font_size=6)
       return nx.draw_networkx(graph, pos=positions, node_color=color_map, hide_ticks=True, node_size=60,font_size=6, with_labels=with_labels, edge_color='#aaaaaa')


   ani = matplotlib.animation.FuncAnimation(fig, update, frames=record.get_size(), interval=interval, repeat=repeat)
   if storage is not None:
       ani.save(filename=storage["filename"], writer=storage["writer"])
   plt.show()


