from flask_wtf import FlaskForm
from flask import session
from wtforms import StringField, SubmitField, SelectField
from wtforms.validators import DataRequired, Length, Email, EqualTo

class MIMPhenoType(FlaskForm):
    MIMPhenoType = StringField('Insert the phenotype MIM number of the disease of interest (e.g. 143100 for Huntington Disease, 104300 for Alzheimer’s Disease, 168600 for Parkinson Disease)', validators =[DataRequired(), Length(min=5, max=20)])
    submit = SubmitField('Create graph and select symptom') 



class Symptoms(FlaskForm):
    symptoms = SelectField(u'Choose a symptom')
