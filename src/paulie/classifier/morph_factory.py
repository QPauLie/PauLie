from paulie.helpers.printing import Debug
from paulie.common.pauli import get_pauli_string, get_pauli_array, multi_pauli_arrays, is_commutate
from paulie.helpers.recording import recording_graph
from paulie.classifier.classification import Morph
import sys, os, traceback

class AppendedException(Exception):
    pass

class DebugException(Exception):
    pass

class DependentException(Exception):
    pass

class NotConnectedException(Exception):
    pass

class MorphFactory(Debug):
      def __init__(self, debug = False, record=None):
          super().__init__(debug)
          self.legs = [] # center is zero leg
          self.record = record
          self.lighting = None
          self.delayed_vertices = []
          self.debug_lighting = None
          self.debug_break = False


      def set_debug(self, debug):
          self.debug = debug

      def recording(self, lighting = None, vertices = None):
           if vertices is None:
               vertices = self.get_vertices()
           if lighting is not None:
              vertices.append(lighting)
           recording_graph(self.record, vertices)

      def set_lighting(self, lighting):
           self.lighting = lighting

      def get_lighting(self):
           return self.lighting

      def get_morph(self):
          return Morph(self.legs)

      def lit(self, lighting, vertix):
          lighting = multi_pauli_arrays(lighting, vertix)
          if self.is_included(lighting):
              raise DependentException()
          return lighting

      def get_lits(self, lighting, vertices=None):
          """Return highlighted vertices (connected to the selected vertex)."""
          if vertices is None:
              vertices = self.get_vertices()
          lits = []
          for v in vertices:
              if v != lighting:
                  if is_commutate(lighting, v) is False:
                      lits.append(v)
          return lits

      def is_empty(self):
          return len(self.legs) == 0

      def is_empty_legs(self):
          return len(self.legs) < 3

      def _find_in_leg(self, leg, v):
          try:
             index = leg.index(v)
          except:
             index = -1
          return index


      def find(self, v):
          for i, leg in enumerate(self.legs):
              index = self._find_in_leg(leg, v)
              if index > -1:
                  return i, index
          return -1, -1

      def is_included(self, v):
          leg_index, vertix_index = self.find(v)
          return leg_index > -1

      def get_vertices(self):
          return [v for leg in self.legs for v in leg]

      def get_center(self):
          if self.is_empty():
              return None
          return self.legs[0][0]

      def set_center(self, v):
          if self.is_empty() is False:
              raise Exception("Center is setted")
          self.legs.append([v])

      def get_long_leg(self):
          if self.is_empty_legs():
              raise Exception("No legs")
          return self.legs[len(self.legs) - 1] 

      def get_one_vertix(self):
          if self.is_empty_legs():
              raise Exception("No legs")
          return self.legs[1][0]

      def _gen_one_legs(self):
          if self.is_empty_legs():
              raise Exception("No legs")
          for i in range(1, len(self.legs)):
              if len(self.legs[i]) == 1:
                  yield self.legs[i]
              else:
                  break

      def get_one_vertices(self):
          vertices = []
          for leg in self._gen_one_legs():
              vertices.append(leg[0].copy())
          return vertices

      def get_PQ(self, lighting):
          one_verices = self.get_one_vertices()
          lits = self.get_lits(lighting, one_verices)
          p = None
          q = None
          for v in one_verices:
              if v in lits:
                  p = v
              else:
                  q = v
              if p is not None and q is not None:
                  return multi_pauli_arrays(p, q), p
          return None, None

      def _gen_two_legs(self):
          if self.is_empty_legs():
              raise Exception("No legs")
          for i in range(1, len(self.legs)):
              if len(self.legs[i]) == 2:
                  yield self.legs[i]
              else:
                  if len(self.legs[i]) > 2:
                       break
      def get_two_legs(self):
          legs = []
          for leg in self._gen_two_legs():
              legs.append((leg[0].copy(), leg[1].copy()))
          return legs

      def get_count_two_legs(self):
          count = 0
          for leg in self._gen_two_legs():
              count += 1
          return count

      def is_two_leg(self):
          count_two_legs = self.get_count_two_legs() 
          if count_two_legs == 0:
              return False

          long_leg = self.get_long_leg()
          if len(long_leg) != 2:
              return True
          return count_two_legs > 1

      def _gen_long_legs(self):
          if self.is_empty_legs():
              raise Exception("No legs")
          is_long = False
          for i in range(len(self.legs)-1, 1, -1):
              if len(self.legs) > 2:
                  yield self.legs[i]
              else:
                  break

      def get_long_legs(self):
          legs = []
          for leg in self._gen_long_legs():
              legs.append(leg.copy())
          return legs


      def append(self, v, lit):
          leg_index, vertix_index = self.find(lit)
          if leg_index == -1:
              raise Exception ("No vertix")
          if leg_index == 0:
              self.legs.insert(1, [v])
              return
          if vertix_index !=  len(self.legs[leg_index]) - 1:
              raise Exception ("The vertix is not the last")

          leg = self.legs[leg_index].copy()
          del self.legs[leg_index]
          leg.append(v)
          if len(leg) >= len(self.legs[len(self.legs) - 1]):
              self.legs.append(leg)
              return
          for i in range(len(self.legs) - 1, 0, -1):
              if len(self.legs[i]) <= len(leg):
                  self.legs.insert(i+1, leg)
                  return

          raise Exception ("Can't append")

      def append_to_center(self, lighting):
          center = self.get_center()
          if len(self.legs) == 1:
              self.append(lighting, center)
              return 

          vertices = self.get_vertices()
          lits = self.get_lits(lighting, vertices)
          if len(lits) == 1:
              if center in lits:
                 self.append(lighting, center)
                 return
              else:
                  lighting = self.lit(lighting, lits[0])
                  lighting = self.lit(lighting, center)
                  self.append(lighting, center)
                  return

          if len(lits) == 2:
              lighting = self.lit(lighting, center)
              self.append(lighting, center)
              return
          raise NotConnectedException()

      def remove(self, v):
          leg_index, vertix_index = self.find(v)
          if leg_index == -1:
              raise Exception ("No vertix")
          if leg_index == 0:
              raise Exception ("Can't delete the center")

          leg = [self.legs[leg_index][i] for i in range(0, vertix_index)]
          del self.legs[leg_index]
          if len(leg) == 0:
              return
          if len(leg) == 1:
             self.legs.insert(1, leg)
             return
          if len(leg) >= len(self.legs[len(self.legs) - 1]):
              self.legs.append(leg)
              return
          for i in range(len(self.legs) - 1, 1, -1):
              if len(self.legs[i]) <= len(leg):
                  self.legs.insert(i+1, leg)
                  return

          raise Exception ("Can't remove")

      def replace(self, v, v_new):
          leg_index, vertix_index = self.find(v)
          if leg_index == -1:
              raise Exception ("No vertix")

          self.legs[leg_index][vertix_index] = v_new

      def merge(self, v, lit):
          if v == lit:
              raise Exception ("Merging into oneself")

          leg_index, vertix_index = self.find(v)
          lit_leg_index, lit_vertix_index = self.find(lit)
          if leg_index == -1 or lit_leg_index == -1:
              raise Exception ("No vertix")
          if lit_leg_index != len(self.legs[lit_leg_index]) - 1:
              raise Exception ("The vertix is not the last")

          if leg_index == lit_leg_index:
              raise Exception ("Merging into self leg")

          lit_leg = self.legs[lit_leg_index].copy()
          leg = self.legs[leg_index].copy()
          if leg_index > lit_leg_index:
              del self.legs[leg_index]
              del self.legs[lit_leg_index]
          else:
              del self.legs[lit_leg_index]
              del self.legs[leg_index]
          for i in range(vertix_index, len(leg)):
              lit_leg.append(leg[i])
          del leg[vertix_index::]

          for i in range(len(self.legs) - 1, 1, -1):
              if len(self.legs[i]) <= len(lit_leg):
                  self.legs.insert(i+1, lit_leg)
                  break

          if len(leg) == 0:
              return
          if len(leg) == 1:
             self.legs.insert(1, leg)
             return
          for i in range(len(self.legs) - 1, 1, -1):
              if len(self.legs[i]) <= len(leg):
                  self.legs.insert(i+1, leg)
                  return
          raise Exception ("Can't merge")

      def print_state(self, lighting = None):
           if self.debug is False:
               return
           self.print_title(f"state graph {len(self.get_vertices())}")
           if lighting is None:
               for i, leg in enumerate(self.legs):
                   if i == 0:
                       self.print_title("center")
                       self.print_vertix(self.legs[i][0])
                       continue
                   self.print_vertices(leg)
               return
           self.print_vertix(lighting, "print state") 
           lits = self.get_lits(lighting)
           for i, leg in enumerate(self.legs):
               if i == 0:
                   self.print_title("center")
                   clit = ""
                   if self.legs[i][0] in lits:
                       clit = "*"
                   self.print_vertix(self.legs[i][0], clit)
                   continue
               self.print_lit_vertices(leg, lits)


      def find_center_and_lits(self, vertices):
          """Find a center with maximum connections. And bring back these connections."""
          center = vertices[0]
          center_lits = []
          for v in vertices:
              lits = self.get_lits(v, vertices)
              if len(lits) > len(center_lits):
                  center = v
                  center_lits = lits
          return center, center_lits


      def append_to_sorted_vertices(self, sorted_vertices, vertices):
          """Adding the next node according to the principle of having a connection with the previous one.
    
          If there are several connected ones, then we insert them after the first one. (To reduce the risk of graph
          reassembly).
          """
          for v in vertices:
              if v in sorted_vertices:
                  vertices.remove(v)
                  continue

              lits = self.get_lits(v, sorted_vertices)
              if len(lits) == 0:
                  continue
              if len(lits) > 1:
                  min_index = len(sorted_vertices)
                  for lit in lits:
                      index = sorted_vertices.index(lit)
                      if index < min_index:
                          min_index = index
                  sorted_vertices.insert(min_index + 1, v)
              else:
                  sorted_vertices.append(v)
              vertices.remove(v)
              return


      def sort_vertices(self, vertices):
          """Sorting nodes in order from the center and then by connections."""
          vertices.sort()
          sorted_vertices = []
          center, lits = self.find_center_and_lits(vertices)
          vertices.remove(center)
          sorted_vertices.append(center)

          for lit in lits:
              vertices.remove(lit)
              if lit not in sorted_vertices:
                  sorted_vertices.append(lit)

          while len(vertices) > 0:
              self.append_to_sorted_vertices(sorted_vertices, vertices)

          return sorted_vertices

      def get_lit_indexes(self, vertices, lits):
          indexes = []
          for i, v in enumerate(vertices):
              if v in lits:
                 indexes.append(i)
          return indexes
       ### pipeline

      def _append_three_graph(self):
          lighting = self.get_lighting()
          self.print_vertix(lighting, "Check and build graph with three vertices")

          if self.is_empty():
              self.set_center(lighting)
              raise AppendedException
 
          if self.is_included(lighting):
              raise DependentException

          if self.is_empty_legs():
              self.append_to_center(lighting)
              raise AppendedException
          self.print_state(lighting)
          self.set_lighting(lighting)
          return self

      def _append_one_legs_in_different_state(self):

          lighting = self.get_lighting()
          self.print_vertix(lighting, "Build graph with different one leg state")
          pq, p = self.get_PQ(lighting)
          if pq is not None:
              lits = self.get_lits(lighting)
              for lit in lits:
                  if lit != p:
                      v = self.lit(pq, lit)
                      self.replace(lit, v)
              self.append(lighting, p)
              long_leg = self.get_long_leg()

              if len(long_leg) > 4:
                  for i in range(4, len(long_leg)):
                      self.append_delayed(long_leg[i])
                  self.remove(long_leg[4])

              raise AppendedException

          self.print_state(lighting)
          self.set_lighting(lighting)
          return self

      def _append_fast(self):
          lighting = self.get_lighting()
          self.print_vertix(lighting, "Append fast")
          center = self.get_center()

          two_legs = self.get_two_legs()
          long_leg = self.get_long_leg()
          if len(long_leg) == 2:
              del two_legs[len(two_legs) - 1]
          if len(two_legs) == 0:
              lits = self.get_lits(lighting)
              if len(lits) == 1:
                  if center in lits:
                      self.append(lighting, center)
                      raise AppendedException
                  long_leg = self.get_long_leg()
                  if long_leg[len(long_leg) - 1] in lits:
                      self.append(lighting, long_leg[len(long_leg) - 1])
                      raise AppendedException

          self.print_state(lighting)
          self.set_lighting(lighting)
          return self

      def _lit_only_long_leg(self):
          lighting = self.get_lighting()

          self.print_vertix(lighting, "Lit only long leg")
          omega = self.get_one_vertix()
          center = self.get_center()
          center_lits = self.get_lits(lighting, [center])
          lits = self.get_lits(lighting, [omega])

          if omega in lits:
              if center not in center_lits:
                  lighting = self.lit(lighting, omega)
              lighting = self.lit(lighting, center)

          two_legs = self.get_two_legs()
          long_leg = self.get_long_leg()

          if len(long_leg) == 2:
              del two_legs[len(two_legs) - 1]

          elif len(long_leg) == 1:
              self.print_state(lighting)
              self.set_lighting(lighting)
              return self
           
          if len(two_legs) == 0:
              self.print_state(lighting)
              self.set_lighting(lighting)
              return self

          long_lits = self.get_lits(lighting, long_leg)

          if len(long_lits) == 0:
              # find lited two leg
              if center in center_lits:
                  lighting = self.lit(lighting, center)
                  lighting = self.lit(lighting, long_leg[0])
                  lighting = self.lit(lighting, omega)
                  lighting = self.lit(lighting, center)
              else:
                  for leg in two_legs:
                      lits = self.get_lits(lighting, leg)
                      v0 = leg[0]
                      v1 = leg[1]
                      if v1 in lits and v0 not in lits:
                          lighting = self.lit(lighting, v1)
                          lits.append(v0)
                      if v0 in lits:
                          lighting = self.lit(lighting, center)
                          lighting = self.lit(lighting, long_leg[0])
                          lighting = self.lit(lighting, omega)
                          lighting = self.lit(lighting, center)
                          break


          # lit second vertix on long leg
          long_lits = self.get_lits(lighting, long_leg)
          lit_indexes = self.get_lit_indexes(long_leg, long_lits)
          
          if 1 not in lit_indexes:
              if 0 in lit_indexes:
                  lighting = self.lit(lighting, long_leg[0])
              else:
                  if len(lit_indexes) == 0:
                      raise NotConnectedException()
                  first_lit = lit_indexes[0]
                  for i in range(first_lit, 1, -1):
                      lighting = self.lit(lighting, long_leg[0])


          long_v0 = long_leg[0]
          long_v1 = long_leg[1]

          for i,leg in enumerate(two_legs):
              lits = self.get_lits(lighting, leg)
              v0 = leg[0]
              v1 = leg[1]
              if v0 in lits and v1 not in lits:
                  lighting = self.lit(lighting, v0)
                  lits.append(v1)
              if v0 not in lits and v1 in lits:
                  lighting = self.lit(lighting, v1)
                  lits.append(v0)

              if v0 in lits and v1 in lits:
                  if center in center_lits:
                      lighting = self.lit(lighting, center)
                      #omega is lited
                      lighting = self.lit(lighting, v1)
                      lighting = self.lit(lighting, v0)
                      lighting = self.lit(lighting, omega)
                      lighting = self.lit(lighting, center)
                  else:
                      long_lits = self.get_lits(lighting, [long_leg[0]])
                      if len(long_lits) == 0:
                          lighting = self.lit(lighting, long_v1)
                      lighting = self.lit(lighting, long_v0)
                      lighting = self.lit(lighting, center)
                      lighting = self.lit(lighting, omega)
                      lighting = self.lit(lighting, v1)
                      lighting = self.lit(lighting, v0)
                      lighting = self.lit(lighting, center)

          self.print_state(lighting)
          self.set_lighting(lighting)
          return self

      def _lit_center(self):
          lighting = self.get_lighting()
          self.print_vertix(lighting, "Liting center")
          center = self.get_center()
          center_lits = self.get_lits(lighting, [center])
          if center not in center_lits:
              long_leg = self.get_long_leg()
              long_lits = self.get_lits(lighting, long_leg)
              lit_indexes = self.get_lit_indexes(long_leg, long_lits)
              first_lit = lit_indexes[0]
              for i in range(first_lit, -1, -1):
                  lighting = self.lit(lighting, long_leg[i])
          self.print_state(lighting)
          self.set_lighting(lighting)
          return self


      def _reduce_long_leg_more_than_one_lits(self):
          lighting = self.get_lighting()
          self.print_vertix(lighting, "Reduce long leg lits")
          omega = self.get_one_vertix()
          center = self.get_center()
          long_leg = self.get_long_leg()
          v0 = long_leg[0]

          while True:
              lits = self.get_lits(lighting, long_leg)
              if len(lits) == 0:
                  self.append(lighting, center)
                  raise AppendedException()




              if len(lits) == 1:
                  if long_leg[0] == lits[0] or long_leg[len(long_leg) - 1] == lits[0]:
                      break

                  if long_leg[0] != lits[0]:
                      lit_indexes = self.get_lit_indexes(long_leg, lits)
                      if lit_indexes[0] < len(long_leg) - 1:
                          for i in range(lit_indexes[0] + 1, len(long_leg)):
                              self.append_delayed(long_leg[i])
                          self.remove(long_leg[lit_indexes[0] + 1])
                      break

              if len(lits) == 2:
                  lit_indexes = self.get_lit_indexes(long_leg, lits)

                  if lit_indexes[0] == 0 and lit_indexes[1] == len(long_leg) - 1:
                      break


              lit_indexes = self.get_lit_indexes(long_leg, lits)
              first = lit_indexes[0]
              second = lit_indexes[1]
              if first > 0 and first + 1 != second:
                  for i in range(second, first - 1, -1):
                      lighting = self.lit(lighting, long_leg[i])
              lits = self.get_lits(lighting, long_leg)
              if len(lits) == 1:
                  if long_leg[0] == lits[0] or long_leg[len(long_leg) - 1] == lits[0]:
                      break

                  if long_leg[0] != lits[0]:
                      lit_indexes = self.get_lit_indexes(long_leg, lits)
                      if lit_indexes[0] < len(long_leg) - 1:
                          for i in range(lit_indexes[0] + 1, len(long_leg)):
                              self.append_delayed(long_leg[i])
                          self.remove(long_leg[lit_indexes[0] + 1])
                      break



              lit_indexes = self.get_lit_indexes(long_leg, lits)
              if len(lit_indexes) > 0:
                  first = lit_indexes[0]
                  if first != 0:
                     for i in range(first, 0, -1):
                         lighting = self.lit(lighting, long_leg[i])

                     lits = self.get_lits(lighting, long_leg)


              lits = self.get_lits(lighting, long_leg)
              if len(lits) == 1:
                  if long_leg[0] == lits[0] or long_leg[len(long_leg) - 1] == lits[0]:
                      break

              if v0 in lits:
                  lighting = self.lit(lighting, center)
                  lighting = self.lit(lighting, omega)
                  lits = self.get_lits(lighting, long_leg)
                  lit_indexes = self.get_lit_indexes(long_leg, lits)
                  first = lit_indexes[0]
                  for i in range(first, -1, -1):
                      lighting = self.lit(lighting, long_leg[i])

                  lighting = self.lit(lighting, center)


          self.print_state(lighting)
          self.set_lighting(lighting)
          return self






      def _append_long_leg_first_and_center_lit(self):

          lighting = self.get_lighting()
          self.print_vertix(lighting, "Append long leg with first lit and center")
          omega = self.get_one_vertix()
          center = self.get_center()
          lits = self.get_lits(lighting, [center, omega])
          is_center_lit = center in lits

          # lit only long leg
          long_leg = self.get_long_leg()
          lits = self.get_lits(lighting, long_leg)

          if is_center_lit and len(lits) == 0:
              #lit only center
              self.print_state(lighting)
              self.append(lighting, center)
              raise AppendedException

          lit_indexes = self.get_lit_indexes(long_leg, lits)

#          self.debugbreak(lighting = lighting, append=True)

          # only long leg and center are lited 
          if len(lit_indexes) == 1 and 0 in lit_indexes:
              # if leg less than 4 or not legs long 2, we can connect to the end of long leg
              is_can_connect_to_end = True
              if self.is_two_leg() and len(long_leg) > 3:
                  is_can_connect_to_end = False

              if is_can_connect_to_end:
                  for v in long_leg:
                      lighting = self.lit(lighting, v)
                  self.append(lighting, long_leg[len(long_leg)-1])
                  raise AppendedException
              else:
                  two_legs = self.get_two_legs()
                  two_leg = two_legs[0]
                  v0 = two_leg[0]
                  v1 = two_leg[1]
                  lighting = self.lit(lighting, center)
                  lighting = self.lit(lighting, v0)
                  lighting = self.lit(lighting, omega)
                  lighting = self.lit(lighting, center)
                  lighting = self.lit(lighting, long_leg[0])
                  lighting = self.lit(lighting, v1)
                  lighting = self.lit(lighting, v0)
                  lighting = self.lit(lighting, center)
                  lighting = self.lit(lighting, long_leg[1])
                  lighting = self.lit(lighting, long_leg[0])
                  lighting = self.lit(lighting, long_leg[2])
                  lighting = self.lit(lighting, long_leg[1])
                  lighting = self.lit(lighting, long_leg[3])
                  lighting = self.lit(lighting, long_leg[2])
                  lighting = self.lit(lighting, omega)
                  lighting = self.lit(lighting, center)
                  lighting = self.lit(lighting, long_leg[0])
                  lighting = self.lit(lighting, long_leg[1])
                  lighting = self.lit(lighting, v0)
                  lighting = self.lit(lighting, v1)
                  lighting = self.lit(lighting, center)
                  lighting = self.lit(lighting, long_leg[0])
                  lighting = self.lit(lighting, v0)
                  lighting = self.lit(lighting, center)
                  self.append(lighting, center)
                  raise AppendedException

          self.print_state(lighting)
          self.set_lighting(lighting)
          return self



      def _append_long_leg_only_last_lit(self):
          lighting = self.get_lighting()
          self.print_vertix(lighting, "Append if long leg last and center are lited")
          center = self.get_center()
          long_leg = self.get_long_leg()

          lits = self.get_lits(lighting, long_leg)
          if len(lits) == 1:

              last_v = long_leg[len(long_leg) - 1]
              if len(long_leg) == 1:
                  lighting = self.lit(lighting, last_v)
                  self.append(lighting, last_v)
                  raise AppendedException
              g = long_leg[len(long_leg) - 2]
              omega = self.get_one_vertix()
              pq = multi_pauli_arrays(omega, lighting)
              new_g = multi_pauli_arrays(pq, g)
              self.remove(last_v)
              self.append(lighting, center)
              self.replace(g, new_g)
              self.append(last_v, lighting)
              long_leg = self.get_long_leg()
              if len(long_leg) > 4:
                  for i in range(4, len(long_leg)):
                      self.append_delayed(long_leg[i])
                  self.remove(long_leg[4])

              raise AppendedException



          self.print_state(lighting)
          self.set_lighting(lighting)
          return self

      def _append_long_leg_last_and_first_lit(self):
          lighting = self.get_lighting()
          self.print_vertix(lighting, "Append if long leg last, first and center are lited")
          omega = self.get_one_vertix()
          center = self.get_center()
          long_leg = self.get_long_leg()
          first_v = long_leg[0]
          for i in range(len(long_leg)-1, 0, -1):
              lighting = self.lit(lighting, long_leg[i])
          lighting = self.lit(lighting, center)
          lighting = self.lit(lighting, omega)
          lighting = self.lit(lighting, first_v)
          lighting = self.lit(lighting, center)
          self.append(lighting, center)
          raise AppendedException


      def _pipeline(self, lighting):
          self.print_vertix(lighting, "Appending vertix to graph")

          # pipeline building
          self.set_lighting(lighting)
          self._append_three_graph()
          self._append_one_legs_in_different_state()
          self._append_fast()
          self._lit_only_long_leg()
          self._lit_center()
          self._reduce_long_leg_more_than_one_lits()
          self._append_long_leg_first_and_center_lit()
          self._append_long_leg_only_last_lit()
          self._append_long_leg_last_and_first_lit()



      def append_delayed(self, v):
           self.delayed_vertices.append(v)

      def restore_delayed(self, vertices):
          for i in range(len(self.delayed_vertices) - 1, -1, -1):
              vertices.insert(0, self.delayed_vertices[i])
          self.delayed_vertices = []
          return vertices

      def set_debug_vertix(self, lighting):
          self.debug_lighting = get_pauli_array(lighting)

      def set_debug_break(self, lighting):
           if lighting == self.debug_lighting:
               self.debug_break = True

      def debugbreak(self, number=None, lighting=None, append=True):
          if self.is_break():
              if append:
                  self.append(lighting, self.get_center())
              self.print_state(lighting)
              raise DebugException()

          if number is not None:
              vertices = self.get_vertices()
              if number <= len(self.get_vertices()):
                  self.print_state(lighting)
                  raise DebugException()
          if lighting is not None:
              self.set_debug_break(lighting)

      def is_break(self):
           return self.debug_break
      #def unlit_two_legs(morph, lighting)
      def build(self, vertices):
          """Transform a connected graph to a cononic type."""
          if len(vertices) == 0:
              return self

          vertices = self.sort_vertices(vertices)
          self.print_vertices(vertices, "init")
          self.recording(vertices=vertices)
          unappended = []
          dependents = []
          #self.set_debug_vertix("ZYIIII")
          while len(vertices) > 0:
              lighting = vertices[0]
              vertices.remove(lighting)
              try:
                  #self.debugbreak(number=13, lighting = lighting)
                  self._pipeline(lighting)
              except AppendedException:
                  vertices = self.restore_delayed(vertices)
                  if lighting in unappended:
                      unappended.remove(lighting)
                  pass
              except DependentException:
                  dependents.append(lighting)
              except NotConnectedException:
                  exc_type, exc_obj, exc_tb = sys.exc_info()
                  self.print_vertix(lighting, f"Exception {traceback.format_exc()}")
                  if lighting not in unappended:
                      unappended.append(lighting)
                      vertices.append(lighting)
              except DebugException:
                  exc_type, exc_obj, exc_tb = sys.exc_info()
                  self.print_vertix(lighting, f"Debug exception {traceback.format_exc()}")
                  break
              except Exception as e:
                  if self.debug:
                      exc_type, exc_obj, exc_tb = sys.exc_info()
                      self.print_vertix(lighting, f"Exception {traceback.format_exc()}")
                  else:
                      self.debuging()
                      self.print_vertix(lighting, f"Exception {e}")
                      self.restore()
                  unappended.append(lighting)


          self.print_state()
          self.print_vertices(unappended, f"unappended")
          self.print_vertices(dependents, f"dependents")
          return self

