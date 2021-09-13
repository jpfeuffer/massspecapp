from pyopenms import FeatureMap

def say(self, msg):
    print ('%s says %s' % (self.name, msg))

if __name__ == '__main__':
    f = FeatureMap()
    FeatureMap.say = say
    f.say('Hello')