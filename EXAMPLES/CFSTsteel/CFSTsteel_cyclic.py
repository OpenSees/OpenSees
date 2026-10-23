import csv
from pathlib import Path

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

def run(output='CFSTsteel_cyclic_python.csv'):
    ops.wipe()

    ops.uniaxialMaterial('CFSTsteel', 1,
        206000.0, 305.0, 0.02, -342.07, 26.58, -74.87, 3.16,
        4.07, 6.38, 3.12, 2.36, 5.82, 1.0, 1.0, 1.0, 1.0, 1.0,
        116.39, 11.81, 1, 4467.68, 38.39)
    ops.testUniaxialMaterial(1)
    last = 0.0
    with Path(output).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['strain', 'stress_MPa', 'tangent_MPa'])
        for peak in (0.0002, -0.0002, 0.003, -0.004, 0.01, -0.02,
                     0.015, -0.03, 0.02, -0.04, 0.01, 0.0):
            for j in range(1, 101):
                strain = last + (peak-last)*j/100.0
                ops.setStrain(strain)
                writer.writerow([strain, ops.getStress(), ops.getTangent()])
            last = peak
    ops.wipe()
    return Path(output)

if __name__ == '__main__':
    print(run())
