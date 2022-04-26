import tkinter as tk

try:
    from pcn.tools.pcn_gui_class import PCNMinerGUI #installed with pip
except:
    try: 
        from pcn_gui_class import PCNMinerGUI #git cloned
    except:
        raise ImportError("PCN-Miner is not correctly installed.")

def main():
    
    window = tk.Tk()
    PCNMinerGUI(window)
    window.mainloop()
    
if __name__ == '__main__':
    main()