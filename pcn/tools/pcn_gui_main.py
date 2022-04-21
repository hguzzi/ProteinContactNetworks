import tkinter as tk
from pcn.tools.pcn_gui_class import PCNMinerGUI

def main():
    
    window = tk.Tk()
    PCNMinerGUI(window)
    window.mainloop()
    
if __name__ == '__main__':
    main()