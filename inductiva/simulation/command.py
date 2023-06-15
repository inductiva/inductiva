"""Base classe for construction of simulator commands.""" 

class Command():
    def __init__(self,
             **kwargs):
        """Command constructor."""
        self.kwargs = kwargs
    
    def get_dict(self):
        """Return the command as a dictionary."""
        return self.kwargs
    

    
