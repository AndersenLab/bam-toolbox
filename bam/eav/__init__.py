from datetime import datetime
from collections import OrderedDict

class eav:
    """
        Entity-Attributes-Value Model
    """
    header_out = False
    def __init__(self, entity, attributes, value, tsv = False, header = False):
        self.entity = entity
        if type(attributes) in [dict, OrderedDict]:
            attributes = [k + "=" + str(v) for k,v in attributes.items()]
        elif type(attributes) != list:
            attributes = [attributes]
        self.attributes = attributes
        self.value = value
        self.tsv = tsv
        if not eav.header_out and header:
            if self.tsv:
                print "ENTITY\t" + '\t'.join([str(x.split("=")[0]) for x in attributes]) + "\tVALUE\tDATE"
                eav.header_out = True
            else:
                print "ENTITY\tATTR\tVALUE\tDATE"
                eav.header_out = True


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
        if self.tsv:
            attributes = '\t'.join(map(str,[x.split("=")[1] for x in self.attributes]))
            out = [self.entity, attributes, self.value, datetime.now()]
            output = map(str, out)
            return '\t'.join(output)
        else:
            # EAV output
            attributes = ';'.join(map(str,self.attributes))
            out = [self.entity, attributes, self.value, datetime.now()]
            output = map(str, out)
            return '\t'.join(output)

def dict_to_eav(entity, ATTR, attributes = {}, tsv = False, header = False):
    for k, value in ATTR.items():
        attr_in = OrderedDict([("ATTR", k)] + attributes.items())
        e = eav(entity, attr_in, value, tsv, header)
        print e

