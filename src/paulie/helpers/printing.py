from paulie.common.pauli import get_pauli_string


def print_node(debug, node, title=""):
    if debug:
        print(f"{title} {get_pauli_string(node)}")

### print array nodes
def print_nodes(debug, nodes, title = ""):
    if debug is False:
        return

    print(f"----{title}----{len(nodes)}---")
    for node in nodes:
        print_node(debug, node)
    print("-------------------")


class Debug:
      def __init__(self, debug):
          self.debug = debug
     
      def print_node(self, node, title=""):
          if self.debug:
              if title != "":
                  print(f"{title}")
              print(f"{get_pauli_string(node)}")
      def print_nodes(self, nodes, title=""):
          if self.debug:
              if title != "":
                  print(f"{title}")
              else:
                  print("--------------------")
              for node in nodes:
                  self.print_node(node)
              print("--------------------")

      def print_title(self, title):
          if self.debug:
              if title != "":
                  print(f"{title}")
