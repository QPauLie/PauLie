import enum

class TypeGraph(enum.Enum):
      A = 0
      B1 = 1
      B2 = 2
      B3 = 3
      NONE = 4

class TypeAlgebra(enum.Enum):
      U = 0
      SU = 1
      SP = 2
      SO = 3


class Morph:
      def __init__(self, legs):
          self.legs = legs # center is zero leg

      def is_empty(self):
          return len(self.legs) == 0

      def is_empty_legs(self):
          return len(self.legs) == 1

      def get_vertices(self):
          return [v for leg in self.legs for v in leg ]

      def counts(self):
          one_legs = 0
          two_legs = 0
          long_vertices = 0
          for i, leg in enumerate(self.legs):
              if i == 0:
                  continue
              if len(leg) == 1:
                  one_legs += 1

              if len(leg) == 2:
                  two_legs += 1
              if len(leg) > 2:
                  if long_vertices > 0:
                      raise Exception("Graph of non-canonical type") 
                  long_vertices += len(leg)

          if long_vertices == 0 and two_legs == 1:
              two_legs = 0
              long_vertices = 2
          if long_vertices == 0 and one_legs > 1 and two_legs == 0:
             one_legs -= 1
             long_vertices = 1

          return one_legs, two_legs, long_vertices

      def get_properties(self):
          if self.is_empty():
              raise Exception("Graph of non-canonical type") 
          if self.is_empty_legs():
              return  TypeGraph.NONE, 0, 0, 0

          one_legs, two_legs, long_vertices = self.counts()
          if two_legs == 0:
              return  TypeGraph.A, one_legs, two_legs, long_vertices
          if  long_vertices == 0:
              return  TypeGraph.B1, one_legs, two_legs, long_vertices
          if  long_vertices == 3:
              return  TypeGraph.B3, one_legs, two_legs, long_vertices
          if  long_vertices == 4:
              return  TypeGraph.B2, one_legs, two_legs, long_vertices
          raise Exception("Graph of non-canonical type") 

      def get_type(self):
          type_graph, one_legs, two_legs, long_vertices = self.get_properties()
          return type_graph

      def get_algebra_properties(self):
          type_graph, one_legs, two_legs, long_vertices = self.get_properties()
          if type_graph == TypeGraph.NONE:
              return TypeAlgebra.U, 1, 1

          if type_graph == TypeGraph.A:
              return TypeAlgebra.SO, one_legs, long_vertices + 2
          if type_graph == TypeGraph.B1:
              return TypeAlgebra.SP, one_legs, 2**two_legs
          if type_graph == TypeGraph.B2:
              return TypeAlgebra.SO, one_legs, 2**(two_legs + 3)
          if type_graph == TypeGraph.B3:
              return TypeAlgebra.SU, one_legs, 2**(two_legs + 2)

      def check_algebra_properties(self, type_algebra = None, nc = None, size = None):
          _type_algebra, _nc, _size = self.get_algebra_properties()
          if type_algebra is None and nc is None and size is None:
              return True
          if type_algebra is not None:
             if type_algebra != _type_algebra:
                 return False
          if nc is not None and size is not None:
              return size == _size and nc == _nc
          if size is not None:
              return size == _size
          if nc is not None:
              return nc == _nc
          return False



class Classification:
      def __init__(self):
          self.morphs = set()

      def add(self, morph):
          self.morphs.add(morph)

      def get_morphs(self):
          return self.morphs

      def get_algebra(self):
          algebras = {}
          for morph in self.morphs:
             type_algebra, nc, size = morph.get_algebra_properties()
             algebra = ""
             if type_algebra == TypeAlgebra.U:
                 algebra = f"u({size})"
             if type_algebra == TypeAlgebra.SO:
                 algebra = f"so({size})"
             if type_algebra == TypeAlgebra.SP:
                 algebra = f"sp({size})"
             if type_algebra == TypeAlgebra.SU:
                 algebra = f"su({size})"
             if algebra in algebras:
                 algebras[algebra] += nc if nc == 1 else 2**(nc-1)
             else:
                 algebras[algebra] = nc if nc == 1 else 2**(nc-1)
          return "+".join([key if v == 1 else str(v) + "*" + key for key, v in algebras.items()])

      def contains_algebra(self, algebra):
          _algebra = self.get_algebra()
          algebra.replace(" ", "")
          return _algebra.find(algebra) > -1

      def is_algebra(self, algebra):
          _algebra = self.get_algebra()
          algebra.replace(" ", "")
          algebras = algebra.split("+")
          _algebras = _algebra.split("+")
          if len(algebras) != len(_algebras):
              return False

          algebras.sort()
          _algebras.sort()
          if i in range(0, len(algebras)):
              if algebras[i] != _algebras[i]:
                  return False
          return True

      def get_subalgebras(self, algebra=None):
          if algebra is None:
              algebra = self.get_algebra()
          else:
              algebra.replace(" ", "")

          return algebra.split("+")


      def get_vertices(self):
          return [v for morph in self.morphs for v in morph.get_vertices() ]



