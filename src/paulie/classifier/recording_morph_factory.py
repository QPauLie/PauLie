from paulie.helpers.printing import Debug
from paulie.helpers.recording import recording_graph
from paulie.classifier.classification import Morph, Classification
import sys, os, traceback

class AppendedException(Exception):
    pass

class DebugException(Exception):
    pass

class DependentException(Exception):
    pass

class NotConnectedException(Exception):
    pass

class RecordingMorphFactory(Debug):
      def __init__(self, debug = False, record=None):
          super().__init__(debug)
          self.legs = [] # center is zero leg
          self.record = record
          self.lighting = None
          self.delayed_vertices = []
          self.debug_lighting = None
          self.debug_break = False
          self.dependents = []


      def set_debug(self, debug):
          self.debug = debug

      #def recording(self, lighting = None, vertices = None):
      #     if vertices is None:
      #         vertices = self.get_vertices()
      #     if lighting is not None:
      #        vertices.append(lighting)
      #     recording_graph(self.record, vertices)

      def set_lighting(self, lighting):
           self.lighting = lighting

      def get_lighting(self):
           return self.lighting

      def get_morph(self):
          return Morph(self.legs, self.dependents)

      def lit(self, lighting, vertix):
          lighting = lighting@vertix
          if self.is_included(lighting):
              recording_graph(self.record, lighting=lighting, dependent=lighting, title=f"Dependent: {lighting}")

              raise DependentException()
          return lighting

      def get_lits(self, lighting, vertices=None):
          """Return highlighted vertices (connected to the selected vertex)."""
          if vertices is None:
              vertices = self.get_vertices()
          return [v for v in vertices if v != lighting and not lighting|v]

#          lits = []
#          for v in vertices:
#              if v != lighting:
#                  if not lighting.commutes_with(v):
#                      lits.append(v)
#          return lits

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
                  return p@q, p, q
          return None, None, None

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
              recording_graph(self.record, lighting=lighting, lits=[center], title=f"Step I: {lighting}")
              recording_graph(self.record, lighting=lighting, lits=[center], appending=center, title=f"Step I: {lighting}")
              self.append(lighting, center)
              return 

          vertices = self.get_vertices()
          lits = self.get_lits(lighting, vertices)
          recording_graph(self.record, lighting=lighting, lits=lits, title=f"Step I: {lighting}")

          if len(lits) == 1:
              if center in lits:
                 self.append(lighting, center)
                 recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), appending=center, title=f"Step I: {lighting}")
                 return                                                                        
              else:
                  lighting = self.lit(lighting, lits[0])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=lits[0], title=f"Step I: {lighting}")
                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step I: {lighting}")
                  self.append(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), appending=center, title=f"Step I: {lighting}")
                  return

          if len(lits) == 2:
              lighting = self.lit(lighting, center)
              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step I: {lighting}")
              self.append(lighting, center)
              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), appending=center, title=f"Step I: {lighting}")
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
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), title=f"Step I: {lighting}")

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
          pq, p, q = self.get_PQ(lighting)
          #recording_graph(self.record, lighting=lighting, lits=[center], appending=center, title=f"StepII: {lighting}")
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), title=f"Step II: {lighting}")

          if pq is not None:
              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), p=p, q=q, title=f"Step II: {lighting}")
              lits = self.get_lits(lighting)
              replacing = []
              for lit in lits:
                  if lit != p:
                      replacing.append(lit)
              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), replacing_vertices=replacing, title=f"Step II: {lighting}")

              for lit in lits:
                  if lit != p:
                      v = pq@lit
                      self.replace(lit, v)
                      replacing.append(lit)
              recording_graph(self.record, collection=self.get_vertices(),lighting=lighting, lits=self.get_lits(lighting), title=f"Step II: {lighting}")

              self.append(lighting, p)
              long_leg = self.get_long_leg()

              if len(long_leg) > 4:
                  removing = []
                  for i in range(4, len(long_leg)):
                      self.append_delayed(long_leg[i])
                      removing.append(long_leg[i]) 
                  if len(removing) > 0:
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), removing_vertices=removing, title=f"Step II: {lighting}")
                  self.remove(long_leg[4])

              recording_graph(self.record, collection=self.get_vertices(), lighting=lighting, lits=self.get_lits(lighting), appending=p, title=f"Step II: {lighting}")

              raise AppendedException

          self.print_state(lighting)
          self.set_lighting(lighting)
          return self

      def _append_fast(self):
          lighting = self.get_lighting()
          self.print_vertix(lighting, "Append fast")
          center = self.get_center()
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), title=f"Step append fast: {lighting}")

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
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), title=f"Step III: {lighting}")

          omega = self.get_one_vertix()
          center = self.get_center()
          center_lits = self.get_lits(lighting, [center])
          lits = self.get_lits(lighting, [omega])

          if omega in lits:
              if center not in center_lits:
                  lighting = self.lit(lighting, omega)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=omega, title=f"Step III: {lighting}")

              lighting = self.lit(lighting, center)
              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step III: {lighting}")

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
              center_lits = self.get_lits(lighting, [center])
              if center in center_lits:
                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step III: {lighting}")

                  lighting = self.lit(lighting, long_leg[0])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[0], title=f"Step III: {lighting}")
                  lighting = self.lit(lighting, omega)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=omega, title=f"Step III: {lighting}")
                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step III: {lighting}")
              else:
                  for leg in two_legs:
                      lits = self.get_lits(lighting, leg)
                      v0 = leg[0]
                      v1 = leg[1]
                      if v1 in lits and v0 not in lits:
                          lighting = self.lit(lighting, v1)
                          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v1, title=f"Step III: {lighting}")
                          lits.append(v0)
                      if v0 in lits:
                          lighting = self.lit(lighting, v0)
                          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v0, title=f"Step III: {lighting}")
                          lighting = self.lit(lighting, center)
                          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step III: {lighting}")
                          lighting = self.lit(lighting, long_leg[0])
                          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[0], title=f"Step III: {lighting}")
                          lighting = self.lit(lighting, omega)
                          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=omega, title=f"Step III: {lighting}")
                          lighting = self.lit(lighting, center)
                          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step III: {lighting}")
                          break


          # lit second vertix on long leg
          long_lits = self.get_lits(lighting, long_leg)
          lit_indexes = self.get_lit_indexes(long_leg, long_lits)
          
          if 1 not in lit_indexes:
              if 0 in lit_indexes:
                  lighting = self.lit(lighting, long_leg[0])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[0], title=f"Step III: {lighting}")
              else:
                  if len(lit_indexes) == 0:
                      raise NotConnectedException()
                  first_lit = lit_indexes[0]
                  for i in range(first_lit, 1, -1):
                      lighting = self.lit(lighting, long_leg[i])
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[i], title=f"Step III: {lighting}")


          long_v0 = long_leg[0]
          long_v1 = long_leg[1]

          for i,leg in enumerate(two_legs):
              lits = self.get_lits(lighting, leg)
              v0 = leg[0]
              v1 = leg[1]
              if v0 in lits and v1 not in lits:
                  lighting = self.lit(lighting, v0)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v0, title=f"Step III: {lighting}")
                  lits.append(v1)
              if v0 not in lits and v1 in lits:
                  lighting = self.lit(lighting, v1)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v1, title=f"Step III: {lighting}")
                  lits.append(v0)

              if v0 in lits and v1 in lits:
                  center_lits = self.get_lits(lighting, [center])
                  if center in center_lits:
                      lighting = self.lit(lighting, center)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step III: {lighting}")
                      #omega is lited
                      lighting = self.lit(lighting, v1)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v1, title=f"Step III: {lighting}")
                      lighting = self.lit(lighting, v0)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v0, title=f"Step III: {lighting}")
                      lighting = self.lit(lighting, omega)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=omega, title=f"Step III: {lighting}")
                      lighting = self.lit(lighting, center)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step III: {lighting}")
                  else:
                      long_lits = self.get_lits(lighting, [long_leg[0]])
                      if len(long_lits) == 0:
                          lighting = self.lit(lighting, long_v1)
                          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_v1, title=f"Step III: {lighting}")
                      lighting = self.lit(lighting, long_v0)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_v0, title=f"Step III: {lighting}")
                      lighting = self.lit(lighting, center)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step III: {lighting}")
                      lighting = self.lit(lighting, omega)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=omega, title=f"Step III: {lighting}")
                      lighting = self.lit(lighting, v1)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v1, title=f"Step III: {lighting}")
                      lighting = self.lit(lighting, v0)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v0, title=f"Step III: {lighting}")
                      lighting = self.lit(lighting, center)
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step III: {lighting}")

          self.print_state(lighting)
          self.set_lighting(lighting)
          return self

      def _lit_center(self):
          lighting = self.get_lighting()
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), title=f"Step IV: {lighting}")
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
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[i], title=f"Step IV: {lighting}")
          self.print_state(lighting)
          self.set_lighting(lighting)
          return self


      def _reduce_long_leg_more_than_one_lits(self):
          lighting = self.get_lighting()
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), title=f"Step IV: {lighting}")
          self.print_vertix(lighting, "Reduce long leg lits")
          omega = self.get_one_vertix()
          center = self.get_center()
          long_leg = self.get_long_leg()
          v0 = long_leg[0]

          while True:
              lits = self.get_lits(lighting, long_leg)
              if len(lits) == 0:
                  self.append(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), appending=center, title=f"Step IV: {lighting}")
                  raise AppendedException()

              if len(lits) == 1:
                  if long_leg[0] == lits[0] or long_leg[len(long_leg) - 1] == lits[0]:
                      break

                  if long_leg[0] != lits[0]:
                      lit_indexes = self.get_lit_indexes(long_leg, lits)
                      if lit_indexes[0] < len(long_leg) - 1:
                          removing = []
                          for i in range(lit_indexes[0] + 1, len(long_leg)):
                              self.append_delayed(long_leg[i])
                              removing.append(long_leg[i])
                          if len(removing) > 0:
                              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), removing_vertices=removing, title=f"Step IV: {lighting}")
                          self.remove(long_leg[lit_indexes[0] + 1])
                          recording_graph(self.record, collection=self.get_vertices(), lighting=lighting, lits=self.get_lits(lighting), title=f"Step IV: {lighting}")
                      break

              if len(lits) == 2:
                  lit_indexes = self.get_lit_indexes(long_leg, lits)

                  if lit_indexes[0] == 0 and lit_indexes[1] == len(long_leg) - 1:
                      break


              lit_indexes = self.get_lit_indexes(long_leg, lits)
              first = lit_indexes[0]
              second = lit_indexes[1]
              if first > 0 and first + 1 != second:
                  for i in range(second, first - 1, -1): ## maybe + 1
                      lighting = self.lit(lighting, long_leg[i])
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[i], title=f"Step IV: {lighting}")

              lits = self.get_lits(lighting, long_leg)
              if len(lits) == 1:
                  if long_leg[0] == lits[0] or long_leg[len(long_leg) - 1] == lits[0]:
                      break

                  if long_leg[0] != lits[0]:
                      lit_indexes = self.get_lit_indexes(long_leg, lits)
                      if lit_indexes[0] < len(long_leg) - 1:
                          removing = []
                          for i in range(lit_indexes[0] + 1, len(long_leg)):
                              self.append_delayed(long_leg[i])
                              removing.append(long_leg[i])
                          if len(removing) > 0:
                              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), removing_vertices=removing, title=f"Step IV: {lighting}")
                          self.remove(long_leg[lit_indexes[0] + 1])
                          recording_graph(self.record, collection=self.get_vertices(), lighting=lighting, lits=self.get_lits(lighting), title=f"Step IV: {lighting}")
                      break



              lit_indexes = self.get_lit_indexes(long_leg, lits)
              if len(lit_indexes) > 0:
                  first = lit_indexes[0]
                  if first != 0:
                     for i in range(first, 0, -1):
                         lighting = self.lit(lighting, long_leg[i])
                         recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[i], title=f"Step IV: {lighting}")

                     lits = self.get_lits(lighting, long_leg)


              lits = self.get_lits(lighting, long_leg)
              if len(lits) == 1:
                  if long_leg[0] == lits[0] or long_leg[len(long_leg) - 1] == lits[0]:
                      break

              if v0 in lits:
                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step IV: {lighting}")
                  lighting = self.lit(lighting, omega)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=omega, title=f"Step IV: {lighting}")
                  lits = self.get_lits(lighting, long_leg)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg, title=f"Step IV: {lighting}")
                  lit_indexes = self.get_lit_indexes(long_leg, lits)
                  first = lit_indexes[0]
                  for i in range(first, -1, -1):
                      lighting = self.lit(lighting, long_leg[i])
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[i], title=f"Step IV: {lighting}")

                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step IV: {lighting}")


          self.print_state(lighting)
          self.set_lighting(lighting)
          return self

      def _append_long_leg_first_and_center_lit(self):

          lighting = self.get_lighting()
          self.print_vertix(lighting, "Append long leg with first lit and center")
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), title=f"Step V: {lighting}")

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
              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), appending=center, title=f"Step V: {lighting}")
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
                      recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v, title=f"Step V: {lighting}")
                  self.append(lighting, long_leg[len(long_leg)-1])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), appending=long_leg[len(long_leg)-1], title=f"Step V: {lighting}")
                  raise AppendedException
              else:
                  two_legs = self.get_two_legs()
                  two_leg = two_legs[0]
                  v0 = two_leg[0]
                  v1 = two_leg[1]
                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, v0)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v0, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, omega)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=omega, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, long_leg[0])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[0], title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, v1)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v1, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, v0)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v0, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, long_leg[1])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[1], title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, long_leg[0])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[0], title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, long_leg[2])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[2], title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, long_leg[1])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[1], title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, long_leg[3])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[3], title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, long_leg[2])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[2], title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, omega)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=omega, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, long_leg[0])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[0], title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, long_leg[1])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[1], title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, v0)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v0, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, v1)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v1, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, long_leg[0])
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[0], title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, v0)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=v0, title=f"Step V: {lighting}")
                  lighting = self.lit(lighting, center)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step V: {lighting}")
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), appending=center, title=f"Step V: {lighting}")
                  self.append(lighting, center)
                  raise AppendedException

          self.print_state(lighting)
          self.set_lighting(lighting)
          return self



      def _append_long_leg_only_last_lit(self):
          lighting = self.get_lighting()
          self.print_vertix(lighting, "Append if long leg last and center are lited")
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), title=f"Step VI: {lighting}")
          center = self.get_center()
          long_leg = self.get_long_leg()

          lits = self.get_lits(lighting, long_leg)
          if len(lits) == 1:

              last_v = long_leg[len(long_leg) - 1]
              if len(long_leg) == 1:
                  lighting = self.lit(lighting, last_v)
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=last_v, title=f"Step VI: {lighting}")
                  recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), appending=last_v, title=f"Step VI: {lighting}")
                  self.append(lighting, last_v)
                  raise AppendedException
              g = long_leg[len(long_leg) - 2]
              omega = self.get_one_vertix()
              pq = omega@lighting
              new_g = pq@g
              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), removing_vertices=[last_v], title=f"Step VI: {lighting}")

              self.remove(last_v)
              recording_graph(self.record, collection=self.get_vertices(), lighting=lighting, lits=self.get_lits(lighting), title=f"Step VI: {lighting}")
              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), appending=center, title=f"Step VI: {lighting}")
              self.append(lighting, center)
              recording_graph(self.record, collection=self.get_vertices(), lighting=lighting, lits=self.get_lits(lighting), title=f"Step VI: {lighting}")
              recording_graph(self.record, lighting=last_v, lits=self.get_lits(lighting), replacing_vertices=[g], title=f"Step VI: {lighting}")
              self.replace(g, new_g)
              recording_graph(self.record, collection=self.get_vertices(), lighting=last_v, lits=self.get_lits(lighting), title=f"Step VI: {lighting}")
              long_leg = self.get_long_leg()
              if len(long_leg) > 4:
                  removing = []
                  for i in range(4, len(long_leg)):
                      self.append_delayed(long_leg[i])
                      removing.append(long_leg[i])
                  if len(removing) > 0:
                      recording_graph(self.record, lighting=last_v, lits=self.get_lits(lighting), removing_vertices=removing, title=f"Step VI: {lighting}")

                  self.remove(long_leg[4])

              recording_graph(self.record, collection=self.get_vertices(), lighting=last_v, lits=self.get_lits(lighting), appending=lighting, title=f"Step VI: {lighting}")
              self.append(last_v, lighting)
              raise AppendedException



          self.print_state(lighting)
          self.set_lighting(lighting)
          return self

      def _append_long_leg_last_and_first_lit(self):
          lighting = self.get_lighting()
          self.print_vertix(lighting, "Append if long leg last, first and center are lited")
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), title=f"Step VII: {lighting}")
          omega = self.get_one_vertix()
          center = self.get_center()
          long_leg = self.get_long_leg()
          first_v = long_leg[0]
          for i in range(len(long_leg)-1, 0, -1):
              lighting = self.lit(lighting, long_leg[i])
              recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=long_leg[i], title=f"Step VII: {lighting}")
          lighting = self.lit(lighting, center)
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step VII: {lighting}")
          lighting = self.lit(lighting, omega)
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=omega, title=f"Step VII: {lighting}")
          lighting = self.lit(lighting, first_v)
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=first_v, title=f"Step VII: {lighting}")
          lighting = self.lit(lighting, center)
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), contracting=center, title=f"Step VII: {lighting}")
          recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting), appending=center, title=f"Step VII: {lighting}")
          self.append(lighting, center)
          raise AppendedException


      def _pipeline(self, lighting):
          self.print_vertix(lighting, "Appending vertix to graph")

          # pipeline building
          self.set_lighting(lighting)
          self._append_three_graph()
          self._append_one_legs_in_different_state()
          #self._append_fast()
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
          self.debug_lighting = lighting

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
      def build(self, generators):
          """Transform a connected graph to a cononic type."""

          if len(generators) == 0:
              return self

          vertices = generators.get_queue().get()
          self.set_debug(generators.get_debug())
          #self.debuging()
          recording_graph(self.record, collection=vertices, title=f"Original graph", init=True)

          self.print_vertices(vertices, "init")
          #self.recording(vertices=vertices)
          unappended = []
          self.dependents = []

          #self.set_debug_vertix("ZYIIII")
          while len(vertices) > 0:
              lighting = vertices[0]
              #if len(self.get_vertices()) > 0:
              recording_graph(self.record, collection=self.get_vertices(), lighting=lighting, title=f"Adding: {lighting}")
              vertices.remove(lighting)
              try:
                  #self.debugbreak(number=4, lighting = lighting)
                  self._pipeline(lighting)
              except AppendedException:
                  vertices = self.restore_delayed(vertices)
                  if lighting in unappended:
                      unappended.remove(lighting)
                  pass
              except DependentException:
                  self.dependents.append(lighting)
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

          classification = Classification()
          classification.add(self.get_morph())

          recording_graph(self.record, collection=self.get_vertices(), title=f"Algebra: {classification.get_algebra()}")

          self.print_state()
          self.print_vertices(unappended, f"unappended")
          self.print_vertices(self.dependents, f"dependents")
          return self

