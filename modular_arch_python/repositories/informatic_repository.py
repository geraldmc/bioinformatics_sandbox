from pydantic import BaseModel
#import uuid

class InformaticsCreate(BaseModel):
  resource: str
  description: str
  has_api: bool = False
  requires_key: bool = False

class Informatics(InformaticsCreate):
  id: int

class InformaticsRepository:
  def __init__(self):
    self.informatics = []
    #self.counter = get_counter()

  def create_informatic(self, informatic: InformaticsCreate) -> Informatics:
    # why does this line appear twice? Once in the repository and once in the service?
    new_informatic = Informatics(id=self.counter.increment(), **informatic.model_dump())
    #new_informatic = Informatics(id=uuid.uuid4().time, **informatic.model_dump())
    self.informatics.append(new_informatic)
    return new_informatic

  def get_informatics(self) -> list[Informatics]:
    return self.informatics

  def get_informatic(self, informatic_id: int) -> Informatics | None:
    for informatic in self.informatics:
      if informatic.id == informatic_id:
        return informatic
    return None

  def update_informatic(self, informatic_id: int, updated_informatic: InformaticsCreate) -> Informatics | None:
    for informatic in self.informatics:
      if informatic.id == informatic_id:
        informatic.title = updated_informatic.title
        return informatic
    return None

  def delete_informatic(self, informatic_id: int) -> bool:
    for index, informatic in enumerate(self.informatics):
      if informatic.id == informatic_id:
        del self.informatics[index]
        return True
    return False
