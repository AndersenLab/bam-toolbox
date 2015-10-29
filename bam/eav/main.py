from datetime import datetime

class eav:
    """
        Entity-Attributes-Value Model
    """
    header = "ENTITY\tATTRIBUTES\tVALUE\tTIMESTAMP"
    def __init__(self, entity, attributes, value):
        self.entity = entity
        if type(attributes) == dict:
            attributes = [k + "=" + str(v) for k,v in attributes.items()]
        elif type(attributes) != list:
            attributes = [attributes]
        self.attributes = attributes
        self.value = value

    def __setattr__(self, name, value):
        # Value is attribute
        if name == "add_attr" or name == "set_attr":
            if type(value) == dict:
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
        attributes = ';'.join(map(str,self.attributes))
        out = [self.entity, attributes, self.value, datetime.now()]
        output = map(str, out)
        return '\t'.join(output)
