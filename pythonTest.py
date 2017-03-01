
from aligner.core import aligner
from aligner.modules import module

a = aligner.Aligner()

a.addModule(module.Printer())

a.step()
