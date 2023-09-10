# ut

import os
import glob
import shutil

SOURCE_PATH = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../'))

def release_on_windows():
    release_source_path = os.path.join(SOURCE_PATH, 'out/build/x64-Clang-Release')
    if not os.path.exists(release_source_path):
        print('Release path not found: ', release_source_path)
        return

    release_dest_path = os.path.join(SOURCE_PATH, 'Release')
    os.makedirs(release_dest_path, exist_ok=True)

    types = ['*.dll', '*.exe', '*.pyd']
    for t in types:
        files = glob.iglob(os.path.join(release_source_path, t))
        for file in files:
            if os.path.isfile(file):
                shutil.copy2(file, release_dest_path)
                print('Copied: ', file)
        print('Release files copied to: ', release_dest_path)

    # run test
    test_data_path = os.path.join(SOURCE_PATH, "tests\data")
    print('Running tests...')
    print("Test data path: ", test_data_path)
    os.system(os.path.join(release_dest_path, 'tests.exe') + f' --data_path "{test_data_path}"')

if __name__ == '__main__':
    if os.name == 'nt':
        release_on_windows()
