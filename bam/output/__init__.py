from datetime import datetime
from collections import OrderedDict

class output_line:
    """
        Entity-Attributes-Value Model
    """
    header_out = False
    def __init__(self, entity, attributes, value, header = False):
        self.entity = entity
        if type(attributes) in [dict, OrderedDict]:
            attributes = [k + "=" + str(v) for k,v in attributes.items()]
        elif type(attributes) != list:
            attributes = [attributes]
        self.attributes = attributes
        self.value = value
        if not output_line.header_out and header:
            print "bam\tcontig\tstart\tend\tproperty\tvalue"
            output_line.header_out = True


    def __setattr__(self, name, value):
        # Value is attribute
        if name == "add_attr" or name == "set_attr":
            if type(value) in [dict, OrderedDict]:
                value = [k + "=" + v for k,v in value.items()]
            elif type(value) != list:
                value = [value]
            if name == "add_attr":
                self.__dict__["attributes"].extend(value)
            else:
                self.__dict__["attributes"] = value
        else:
            self.__dict__[name] = value


    def __repr__(self):
        attributes = '\t'.join(map(str,[x.split("=")[1] for x in self.attributes]))
        out = [self.entity, attributes, self.value]
        output = map(str, out)
        return '\t'.join(output)
       
