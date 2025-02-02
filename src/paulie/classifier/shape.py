
      
from paulie.helpers.printing import Debug
from paulie.common.pauli import get_pauli_string


class Shape(Debug):                                                     
      def __init__(self, debug = False):
          super().__init__(debug)
          self.center = None
          self.prohibitedLine = False
          self.reset()

      def reset(self):
          self.legs = {1:[], 2: [], 0:[]}
          self.lenLine = 0
          self.links = {}
          self.directlinks = {}
          self.prohibited = None
          self.prohibitedsForPop = []

      def get_center(self):
          return self.center

      def set_prohibited_line(self, prohibitedLine):
           self.prohibitedLine = prohibitedLine

      def is_prohibited_line(self):
          return self.prohibitedLine

      def append_prohibiteds_for_pop(self, v):
          self.prohibitedsForPop.append(v)

      def clear_prohibiteds_for_pop(self):
          self.prohibitedsForPop = []

      def is_prohibiteds_for_pop(self, v):
          return v in self.prohibitedsForPop

      def get_len_line(self):
          return self.lenLine

      def set_center(self, center):
          self.center = center

      def print_state(self, title):
          self.print_title(title)
          self.print_node(self.center, "center")
          self.print_nodes(self.legs[1], "one")
          self.print_nodes(self.legs[2], "two")
          self.print_nodes(self.legs[0], f"line {self.lenLine}")

      def reset_prohibited(self):
          self.prohibited = None

      def set_prohibited(self, prohibited):
          self.prohibited = prohibited

      def is_prohibited(self, v):
          return self.prohibited == v

      def get_lit(self, v):
           return self.links[get_pauli_string(v)]

      def append_is_canonical(self, v, lit):
          if self.is_prohibited(lit):
              return False

          if v == self.center:
               self.reset_prohibited()
               return True

          if lit == self.center:
              if self.lenLine > 4:
                  return False
              #if len(self.legs[2]) == 0 and len(self.legs[1]) > 2:
              #    return False
              self.legs[1].append(v)
              self.links[get_pauli_string(v)] = lit
              self.reset_prohibited()
              return True
          if lit in self.legs[1]:
              if len(self.legs[1]) == 1:
                  return False
              self.legs[1].remove(lit)
              self.legs[2].append(v)
              self.links[get_pauli_string(v)] = lit
              self.directlinks[get_pauli_string(lit)] = v
              self.reset_prohibited()
              return True
          if lit in self.legs[2] and len(self.legs[0]) == 0:
              self.legs[2].remove(lit)
              self.legs[0].append(v)
              self.lenLine +=3
              self.links[get_pauli_string(v)] = lit
              self.directlinks[get_pauli_string(lit)] = v
              self.reset_prohibited()
              return True
          if lit in self.legs[0]:
              if len(self.legs[2]) > 0 and self.lenLine == 4:
                  return False
              if self.lenLine == 4 and self.is_prohibited_line():
                  return False
              self.print_title(f"appending to canonic {self.lenLine}")
              self.legs[0].remove(lit)
              self.legs[0].append(v)
              self.links[get_pauli_string(v)] = lit
              self.directlinks[get_pauli_string(lit)] = v
              self.lenLine +=1
              self.reset_prohibited()
              return True
          return False

      def pop_prev_vertix(self, v):
          lit = self.links[get_pauli_string(v)]
          if get_pauli_string(lit) in self.directlinks:
              del self.directlinks[get_pauli_string(lit)]
          del self.links[get_pauli_string(v)]
          return lit

      def find_vertix_for_pop(self):
          if len(self.legs[0]) > 0:
              if self.is_prohibiteds_for_pop(self.legs[0][0]) is False:
                  return 0, self.legs[0][0]
          for v in self.legs[2]:
              if self.is_prohibiteds_for_pop(v) is False:
                  return 2, v
          for v in self.legs[1]:
              if self.is_prohibiteds_for_pop(v) is False:
                  return 1, v
          return -1, None

      def pop_vertix(self):
          leg, v = self.find_vertix_for_pop()
          if leg == -1:
              return None, None
          if leg == 0:
              lit = self.pop_prev_vertix(v)
              self.lenLine -= 1
              self.legs[0].remove(v)
              if self.lenLine == 2:
                  self.legs[2].append(lit)
                  self.lenLine = 0
              else:
                  self.legs[0].append(lit)
                  #self.lenLine += 1
              self.set_prohibited(lit)
              return v, lit
          else:
              if leg == 2:
                  lit = self.pop_prev_vertix(v)
                  self.legs[2].remove(v)
                  self.legs[1].append(lit)
                  self.set_prohibited(lit)
                  return v, lit
              else: 
                if leg == 1:
                    lit = self.pop_prev_vertix(v)
                    self.legs[1].remove(v)
                    self.set_prohibited(lit)
                    return v, lit
                else:
                    return None, None

      def remove_vertix(self, v):
          if v in self.legs[0]:
              self.legs[0].remove(v)
              lit = self.pop_prev_vertix(v)
              self.lenLine -= 1
              if self.lenLine == 2:
                  self.legs[2].append(lit)
                  self.lenLine = 0
              else:
                  self.legs[0].append(lit)
                  #self.lenLine += 1
          if v in self.legs[1]:
              self.legs[1].remove(v)
          if v in self.legs[2]:
              self.legs[2].remove(v)
              lit = self.pop_prev_vertix(v)
              self.legs[1].append(lit)

          if get_pauli_string(v) in self.links:
             lit = self.links[get_pauli_string(v)]
             self.set_prohibited(lit)
             del self.links[get_pauli_string(v)]
             if get_pauli_string(lit) in self.directlinks:
                 del self.directlinks[get_pauli_string(lit)]

      def is_end(self, v):
          if v in self.legs[1]:
              return True
          if v in self.legs[2]:
              return True
          if v in self.legs[0]:
              return True
          return False

      def get_route_to_end(self, v):
          if v == self.center:
              return []
          if self.is_end(v):
              return [v]
          w = v
          route = [w]
          while get_pauli_string(w) in self.directlinks:
               w = self.directlinks[get_pauli_string(w)]
               route.insert(0, w)
          return route

      def decrease_len_line(self):
          vertixs = []
          if self.lenLine == 0 and len(self.legs[1]) == 0:
              return vertixs
          if self.lenLine > 4:
              while(self.lenLine > 4):
                  v = self.legs[0][0]
                  self.remove_vertix(v)
                  vertixs.insert(0, v)
              self.set_prohibited(self.legs[0][0]) 
          else:
              if self.lenLine < 5:
                  if self.lenLine > 0:
                      v = self.legs[0][0]
                  else:
                      if len(self.legs[2]) > 0:
                          v = self.legs[2][0]
                      else:
                          v = self.legs[1][0]
                  lit = self.get_lit(v)
                  self.remove_vertix(v)
                  vertixs.insert(0, v)
                  self.set_prohibited(lit) 
          return vertixs

      def get_type(self):
          if  len(self.legs[2]) == 0 or (len(self.legs[2]) == 1 and self.lenLine == 0):
              lenLine = self.lenLine
              if lenLine == 0 and len(self.legs[2]) == 1:
                 lenLine = 2
              if lenLine == 0 and len(self.legs[1]) > 0:
                 lenLine = 1
              if lenLine > 0:
                  return "A", lenLine + 1 , len(self.legs[1]) - 1, 0
              else:
                  return "None", lenLine, 0, 0

          if self.lenLine == 0 and len(self.legs[2]) > 0:
              return "B1", 0, len(self.legs[1]) - 1, len(self.legs[2])
          if self.lenLine == 3:
              return "B2", 3, len(self.legs[1]) - 1, len(self.legs[2])
          if self.lenLine == 4:
              return "B3", 4, len(self.legs[1]) - 1, len(self.legs[2])
          return "None", 0, 0, 0

      def get_algebra(self):
          type_graph, nl, nc, n2 = self.get_type()
          #print(f"type={type_graph}, nl = {nl}, nc={nc}, n2={n2}")
          n = 2**nc
          algebras = []
          for index in range(0, n):
              if type_graph == "A":
                 algebras.append(f"so({nl + 1})")
              if type_graph == "B1":
                 algebras.append(f"sp({2**n2})")
              if type_graph == "B2":
                 algebras.append(f"so({2**(n2+3)})")
              if type_graph == "B3":
                 algebras.append(f"su({2**(n2+2)})")
              if type_graph == "None":
                 algebras.append("u(1)")

          return " + ".join(algebras)