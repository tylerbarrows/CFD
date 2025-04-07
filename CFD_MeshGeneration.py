def meshGenerationFunction(geometryInputs):
    # TO ADD
    # error when diverging is selected and geometry inputs show converging
    # error when inputs are out of bounds
    # add 2D option eventually
    # be cool

    # import
    import math
    import numpy as np

    class meshGeneration_result:
        # importing variables
        shape = geometryInputs.shape
        dimension = geometryInputs.dimension
        length = geometryInputs.length
        heightInlet = geometryInputs.heightInlet
        heightThroat = geometryInputs.heightThroat
        heightOutlet = geometryInputs.heightOutlet
        dx = geometryInputs.dx
        dy = geometryInputs.dy

        # flag variable
        flag = 0 # all good if flag = 0

        # find x_vector, h_vector, volume_vector
        match dimension:
            case '1D':
                match shape:
                    case 'straight':
                        dx = length/math.ceil(length/dx)
                        x_vector = np.arange(0,length,dx)
                        h_vector = np.linspace(heightInlet,heightInlet,len(x_vector))
                    case 'converging':
                        dx = length/math.ceil(length/dx)
                        x_vector = np.arange(0,length,dx)
                        h_vector = np.linspace(heightInlet,heightOutlet,len(x_vector))
                    case 'diverging':
                        dx = length/math.ceil(length/dx)
                        x_vector = np.arange(0,length,dx)
                        h_vector = np.linspace(heightInlet,heightOutlet,len(x_vector))
                    case 'converging-diverging':
                        dx = length/math.ceil(length/dx)
                        x_vector = np.arange(0,length,dx)
                        h_vector = np.append(np.linspace(heightInlet,heightThroat,math.floor(len(x_vector)/2)),np.linspace(heightThroat+(heightOutlet-heightThroat)/(len(x_vector)-math.floor(len(x_vector)/2))*1,heightOutlet,len(x_vector)-math.floor(len(x_vector)/2)))
                        # h_vector = np.append(np.linspace(heightInlet,heightThroat,math.floor(len(x_vector)/2)),np.linspace(heightThroat,heightOutlet,len(x_vector)-math.floor(len(x_vector)/2)))
                    case _:
                        flag = 1
            case '2D':
                match shape:
                    case 'straight':
                        dx = 1
                    case 'converging':
                        dx = 1
                    case 'diverging':
                        dx = 1
                    case 'converging-diverging':
                        dx = 1
                    case _:
                        flag = 1
            case _:
                flag = 1

    return meshGeneration_result